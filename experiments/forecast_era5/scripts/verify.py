#!/usr/bin/env python3
"""Verify the PlaSiC forecasts against ERA5 and store metrics/maps.

中文说明：用 ERA5 检验 PlaSiC 预报并保存评分指标与场图。

Model output comes from hourly first-day and 6-hourly days 2–10 binary frame
streams; the verification truth is the ERA5 analysis on the same Gaussian grid
and standard pressure levels.  Scores use Gaussian quadrature weights in
latitude and uniform longitude weights, matching the area-weighting principle
used by WeatherBench 2.

中文说明：模式输出来自首日逐小时、第 2–10 天每 6 小时的二进制帧流；检验
真值是插到同一高斯网格和标准气压层上的 ERA5 分析场。评分采用纬度方向的高斯
求积权重和经度方向的均匀权重，与 WeatherBench 2 的面积加权原则一致。

Time matching is strict: the first 24 h are verified against the hourly ERA5
"shock" files, later lead times against the 6-hourly verification files.  A
frame whose verification time has no ERA5 counterpart within 30 minutes is
skipped (metric = NaN).

中文说明：时间匹配非常严格：前 24 小时使用逐小时 ERA5 “shock” 文件检验，
之后的时效使用 6 小时检验文件。如果某帧的检验时刻在 30 分钟内找不到对应的
ERA5 时刻，则跳过该帧（指标记为 NaN）。

Outputs, per case and method, a compressed ``.npz`` holding lead times, RMSE
and bias arrays per variable/level, and the fields at selected map times.

中文说明：对每个个例和初值方法输出一个压缩 ``.npz`` 文件，内含预报时效、
逐变量/层次的 RMSE 与偏差数组，以及选定地图时刻的场。

Use ``--require-exact`` when the six-hour ERA5 archive has been downloaded from
CDS.  In that mode the verifier refuses interpolated or incomplete truth data,
and also refuses stale model frame streams that do not contain every expected
hourly/6-hourly lead.

中文说明：当 6 小时 ERA5 归档已从 CDS 下载完成时使用 ``--require-exact``。
该模式下检验脚本拒绝插值或不完整的真值数据，也拒绝缺少任何预期逐小时/
6 小时时效的过期模式帧流。
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_ic  # noqa: E402
import common  # noqa: E402
import frame_io  # noqa: E402
import restart_io  # noqa: E402

VERIF_ROOT = common.DATA_ROOT / "verification"  # 检验输出目录 / verification output directory
MAP_TIMES_H = [0, 24, 48, 72, 120, 168, 240]    # 保存场图的预报时效 / lead times kept as maps
TIME_TOLERANCE = np.timedelta64(30, "m")        # 宽松模式下的匹配容差 / match tolerance in compatibility mode
ZERO_TOLERANCE = np.timedelta64(0, "m")         # 严格模式下要求精确匹配 / exact match in strict mode

# 中文说明：存档的高层检验变量与层次（T500 等由该表驱动）。
# English: archived upper-air verification variables and their levels.
UPPER_LEVELS = {
    "t": [850, 500, 200],
    "u": [850, 250],
    "v": [850, 250],
    "q": [850, 700],
}

G = common.GRAVITY   # 重力加速度 / gravitational acceleration
RD = common.GASCON   # 干空气气体常数 / gas constant of dry air


class CoverageError(ValueError):
    """Raised when strict verification cannot prove complete time coverage.

    中文说明：当严格检验无法证明时间覆盖完整时抛出。
    """


def expected_times(t0: np.datetime64, lead_hours: "np.ndarray | list[int]") -> np.ndarray:
    """Return exact UTC timestamps for a case's requested lead hours.

    中文说明：返回某个例所要求预报时效对应的精确 UTC 时刻。
    """
    return np.asarray(
        [t0 + np.timedelta64(int(hour), "h") for hour in lead_hours],
        dtype="datetime64[ns]",
    )


def dataset_times(dataset: xr.Dataset) -> np.ndarray:
    """Return a normalized, nanosecond-resolution valid-time coordinate.

    中文说明：返回归一化为纳秒精度的有效时间坐标。
    """
    if "valid_time" not in dataset:
        raise CoverageError("ERA5 dataset has no valid_time coordinate")
    return np.asarray(dataset["valid_time"].values, dtype="datetime64[ns]")


def missing_times(dataset: xr.Dataset, expected: np.ndarray) -> np.ndarray:
    """Return expected timestamps absent from a dataset, using exact equality.

    中文说明：按精确相等判断，返回数据集中缺失的预期时刻。
    """
    # ``datetime64[ns].tolist()`` becomes integer nanoseconds on recent
    # NumPy versions, so comparing those values directly with datetime64
    # scalars would falsely report every timestamp as missing.
    # 中文说明：在新版 NumPy 中 ``datetime64[ns].tolist()`` 会变成整数纳秒，
    # 直接与 datetime64 标量比较会把所有时刻都误报为缺失，因此这里统一转成
    # int64 纳秒后再做集合判断。
    available = dataset_times(dataset).astype("datetime64[ns]").astype("int64")
    requested = np.asarray(expected, dtype="datetime64[ns]")
    requested_int = requested.astype("int64")
    return requested[~np.isin(requested_int, available)]


def require_coverage(dataset: xr.Dataset, expected: np.ndarray, label: str) -> None:
    """Raise a detailed error when any expected timestamp is absent.

    中文说明：只要有一个预期时刻缺失，就抛出带详细信息的错误。
    """
    missing = missing_times(dataset, expected)
    if missing.size:
        available = dataset_times(dataset)
        raise CoverageError(
            f"{label} is incomplete: missing {missing.size} expected timestamps "
            f"({format_times(missing)}); available range is "
            f"{format_times(available[:1])} to {format_times(available[-1:])}."
        )


def format_times(values: np.ndarray, limit: int = 8) -> str:
    """Format timestamps for error messages (at most ``limit`` + summary).

    中文说明：把时刻列表格式化为错误信息（最多显示 ``limit`` 个，其余给出总数）。
    """
    values = np.asarray(values, dtype="datetime64[m]")
    labels = [str(value) for value in values[:limit]]
    if values.size > limit:
        labels.append(f"... ({values.size} total)")
    return ", ".join(labels) if labels else "none"


def open_pl(paths: "list[Path]") -> xr.Dataset:
    """Open pressure-level files, concatenating along ``valid_time``.

    中文说明：打开多个气压层文件，并沿 ``valid_time`` 维拼接、去除重复时刻。
    """
    if not paths:
        raise FileNotFoundError("no ERA5 pressure-level files found")
    datasets = [xr.open_dataset(path) for path in paths]
    datasets = [_normalize_time_coordinate(dataset) for dataset in datasets]
    if len(datasets) == 1:
        return datasets[0]
    ds = xr.concat(datasets, dim="valid_time")
    _, index = np.unique(ds["valid_time"].values, return_index=True)
    return ds.isel(valid_time=np.sort(index))


def open_time_datasets(paths: "list[Path]") -> xr.Dataset:
    """Open daily ERA5 files without a day/time Cartesian product.

    中文说明：打开按天拆分的 ERA5 文件，避免 day×time 笛卡尔积问题。
    """
    paths = sorted(paths)
    if not paths:
        raise FileNotFoundError("no ERA5 verification files found")
    datasets = [_normalize_time_coordinate(xr.open_dataset(path)) for path in paths]
    if len(datasets) == 1:
        return datasets[0]
    dim = "valid_time"
    ds = xr.concat(datasets, dim=dim)
    _, index = np.unique(ds[dim].values, return_index=True)
    return ds.isel({dim: np.sort(index)})


def _normalize_time_coordinate(dataset: xr.Dataset) -> xr.Dataset:
    """Normalize CDS ``time``/``valid_time`` variants to ``valid_time``.

    中文说明：把 CDS 的 ``time``/``valid_time`` 两种写法统一为 ``valid_time``。
    """
    if "valid_time" in dataset:
        return dataset
    if "time" in dataset:
        return dataset.rename({"time": "valid_time"})
    raise CoverageError("ERA5 dataset has neither valid_time nor time")


# ---------------------------------------------------------------------------
# ERA5 truth
# ERA5 真值
# ---------------------------------------------------------------------------
class Truth:
    """Lazy access to all ERA5 truth fields needed by the verification.

    中文说明：按需提供检验所需的全部 ERA5 真值场。

    In strict mode every expected six-hourly/hourly timestamp must exist.  In
    compatibility mode missing timestamps are filled by linear time
    interpolation so that older archives can still be evaluated.
    中文说明：严格模式下所有预期的 6 小时/逐小时时刻都必须存在；兼容模式下
    对缺失时刻做线性时间插值，以便仍能评估较早的归档。
    """

    def __init__(self, case: common.ForecastCase, require_exact: bool = False):
        self.case = case
        self.require_exact = bool(require_exact)
        self.t0 = np.datetime64(common.init_datetime(case))
        directory = common.ERA5_ROOT / case.key
        self.expected_six_hour_times = expected_times(self.t0, np.arange(0, 241, 6))
        self.expected_hourly_times = expected_times(self.t0, np.arange(0, 25, 1))
        # 优先使用新的 6 小时/逐小时归档，旧文件仅作兼容回退。
        # Prefer the new 6-hourly/hourly archives; legacy files are fallbacks.
        verif_pl6 = sorted(directory.glob("verif6_pl*.nc"))
        verif_sfc6 = sorted(directory.glob("surface_verif6*.nc"))
        verif_pl = verif_pl6 or sorted(directory.glob("verif_pl*.nc"))
        verif_sfc = verif_sfc6 or sorted(directory.glob("surface_verif*.nc"))
        self.pl = open_time_datasets(verif_pl)
        self.sfc = open_time_datasets(verif_sfc)
        pressure_missing = missing_times(self.pl, self.expected_six_hour_times)
        surface_missing = missing_times(self.sfc, self.expected_six_hour_times)
        if self.require_exact and pressure_missing.size:
            require_coverage(self.pl, self.expected_six_hour_times, "ERA5 pressure-level verification")
        if self.require_exact and surface_missing.size:
            require_coverage(self.sfc, self.expected_six_hour_times, "ERA5 surface verification")
        self.pressure_interpolated = bool(pressure_missing.size)
        self.surface_interpolated = bool(surface_missing.size)
        # Compatibility mode can still evaluate the archived experiment by
        # constructing virtual six-hour targets. Strict mode never enters
        # this branch and therefore never interpolates truth in time.
        # 中文说明：兼容模式通过构造“虚拟的”6 小时目标，仍可评估已归档的试验；
        # 严格模式不会进入该分支，因此永远不会在时间上插值真值。
        if self.pressure_interpolated and not self.require_exact:
            self.pl = self.pl.interp(valid_time=self.expected_six_hour_times)
        if self.surface_interpolated and not self.require_exact:
            self.sfc = self.sfc.interp(valid_time=self.expected_six_hour_times)
        self.verification_source = (
            "exact six-hourly ERA5 timestamps"
            if not (self.pressure_interpolated or self.surface_interpolated)
            else "linear interpolation of incomplete archived verification timestamps"
        )
        self.pressure_files = tuple(str(path) for path in verif_pl)
        self.surface_files = tuple(str(path) for path in verif_sfc)
        # 初始地表文件同样兼容按天拆分与单文件两种形式。
        # Initial surface files support both the daily split and single-file layouts.
        ic_paths = sorted(directory.glob("surface_ic_*.nc"))
        if not ic_paths:
            ic_paths = sorted(directory.glob("surface_ic.nc"))
        self.ic_sfc = open_time_datasets(ic_paths)
        shock_paths = sorted(directory.glob("shock_pl_*.nc"))
        actual_shock = open_pl(shock_paths) if shock_paths else None
        shock_missing = (
            missing_times(actual_shock, self.expected_hourly_times)
            if actual_shock is not None else self.expected_hourly_times
        )
        if self.require_exact and actual_shock is None:
            raise CoverageError(
                f"{case.key}: hourly pressure-level shock files are required; "
                "none were found."
            )
        if self.require_exact and shock_missing.size:
            require_coverage(actual_shock, self.expected_hourly_times,
                             "ERA5 hourly pressure-level shock")
        if self.require_exact:
            missing_variables = sorted(set(UPPER_LEVELS) - set(actual_shock.data_vars))
            if missing_variables:
                raise CoverageError(
                    "ERA5 hourly pressure-level shock is missing variables: "
                    + ", ".join(missing_variables)
                )
        hourly_fallback = (
            self.pl.interp(valid_time=self.expected_hourly_times)
            if not self.require_exact else None
        )
        if actual_shock is None:
            self.shock = hourly_fallback
            self.shock_fallback = None
        else:
            # Keep the archived hourly stream (which may use 32 pressure
            # levels) and retain a separate fallback for variables it lacks.
            # 中文说明：保留归档的逐小时流（可能使用 32 个气压层），同时为它
            # 缺少的变量保留独立的回退数据。
            self.shock = actual_shock
            self.shock_fallback = hourly_fallback
        # New downloads are split by day to preserve the requested hourly
        # axis.  Retain the single-file path as a compatibility fallback for
        # older archives, but never mix it with the exact daily files.
        # 中文说明：新的下载按天拆分，以保留所请求的逐小时时间轴；单文件路径
        # 仅作为旧归档的兼容回退，且绝不与精确的按天文件混用。
        shock_surface_paths = sorted(directory.glob("shock_surface_*.nc"))
        if shock_surface_paths:
            actual_shock_sfc = open_time_datasets(shock_surface_paths)
        else:
            shock_surface = directory / "shock_surface.nc"
            actual_shock_sfc = (
                _normalize_time_coordinate(xr.open_dataset(shock_surface))
                if shock_surface.exists() else None
            )
        shock_surface_missing = (
            missing_times(actual_shock_sfc, self.expected_hourly_times)
            if actual_shock_sfc is not None else self.expected_hourly_times
        )
        if self.require_exact and actual_shock_sfc is None:
            raise CoverageError(
                f"{case.key}: hourly surface shock file is required; none was found."
            )
        if self.require_exact and shock_surface_missing.size:
            require_coverage(actual_shock_sfc, self.expected_hourly_times,
                             "ERA5 hourly surface shock")
        if self.require_exact:
            missing_variables = sorted({"sp", "t2m"} - set(actual_shock_sfc.data_vars))
            if missing_variables:
                raise CoverageError(
                    "ERA5 hourly surface shock is missing variables: "
                    + ", ".join(missing_variables)
                )
        if actual_shock_sfc is None:
            self.shock_sfc = (
                self.sfc.interp(valid_time=self.expected_hourly_times)
                if not self.require_exact else None
            )
            self.shock_sfc_fallback = None
        else:
            self.shock_sfc = actual_shock_sfc
            self.shock_sfc_fallback = (
                self.sfc.interp(valid_time=self.expected_hourly_times)
                if not self.require_exact else None
            )
        self.shock_pressure_files = tuple(str(path) for path in shock_paths)
        self.shock_surface_file = ",".join(str(path) for path in shock_surface_paths)
        if not self.shock_surface_file:
            legacy_shock_surface = directory / "shock_surface.nc"
            self.shock_surface_file = (
                str(legacy_shock_surface) if legacy_shock_surface.exists() else ""
            )
        # 记录各类缺测情况，便于报告数据覆盖率 / record coverage gaps for reporting
        self.coverage = {
            "pressure_missing": pressure_missing,
            "surface_missing": surface_missing,
            "shock_pressure_missing": shock_missing,
            "shock_surface_missing": shock_surface_missing,
        }
        # 严格模式要求精确匹配；兼容模式允许 ±30 分钟容差。
        # Strict mode demands exact matches; compatibility mode allows +/-30 min.
        self.match_tolerance = ZERO_TOLERANCE if self.require_exact else TIME_TOLERANCE
        self._cache: dict = {}

    def _match(self, dataset: xr.Dataset | None, when: np.datetime64) -> int:
        """Index of the closest time within tolerance, or -1.

        中文说明：返回容差范围内最接近 ``when`` 的时间下标；找不到返回 -1。
        """
        if dataset is None:
            return -1
        times = dataset_times(dataset)
        index = int(np.argmin(np.abs(times - when)))
        return index if abs(times[index] - np.datetime64(when, "ns")) <= self.match_tolerance else -1

    def _dataset(self, when: np.datetime64,
                 name: str | None = None) -> tuple[xr.Dataset | None, int]:
        """Pick the dataset whose time axis matches ``when`` (strict).

        中文说明：选择时间轴与 ``when`` 匹配的数据集（严格匹配）。

        The hourly "shock" archive is preferred inside its window (when it
        carries the requested variable), with a fallback to the 6-hourly
        archive at its edges.

        中文说明：在逐小时 “shock” 归档覆盖的窗口内优先使用它（前提是其中
        含有所需变量）；窗口边缘则回退到 6 小时归档。
        """
        shock_dataset = self.shock
        if (shock_dataset is not None and name is not None
                and name not in shock_dataset.data_vars
                and self.shock_fallback is not None):
            shock_dataset = self.shock_fallback
        in_shock_window = (
            shock_dataset is not None
            and self.t0 <= when <= self.t0 + np.timedelta64(24, "h")
            and (name is None or name in shock_dataset.data_vars)
        )
        if in_shock_window:
            index = self._match(shock_dataset, when)
            if index >= 0:
                return shock_dataset, index
        return self.pl, self._match(self.pl, when)

    def _surface_dataset(self, when: np.datetime64) -> tuple[xr.Dataset | None, int]:
        """Surface counterpart of :meth:`_dataset`.

        中文说明：与 :meth:`_dataset` 对应的地表场版本。
        """
        shock_dataset = self.shock_sfc
        in_shock_window = (
            shock_dataset is not None
            and self.t0 <= when <= self.t0 + np.timedelta64(24, "h")
        )
        if in_shock_window:
            index = self._match(shock_dataset, when)
            if index >= 0:
                return shock_dataset, index
        return self.sfc, self._match(self.sfc, when)

    def upper(self, name: str, level: int, when: np.datetime64):
        """Regridded upper-air truth field (or ``None`` when unavailable).

        中文说明：返回插值到模式网格的高层真值场；无对应数据时返回 ``None``。

        ``name`` uses the ERA5 short names (t/u/v/q/z) and ``level`` is in hPa.
        中文说明：``name`` 使用 ERA5 短名（t/u/v/q/z），``level`` 单位为 hPa。
        """
        key = (name, level, str(when))
        if key in self._cache:
            return self._cache[key]
        dataset, index = self._dataset(when, name)
        if index < 0:
            return None
        selected = dataset.isel(valid_time=index)
        # 选择最近的气压层 / select the nearest pressure level
        level_index = int(np.argmin(np.abs(selected["pressure_level"].values - level)))
        field = selected[name].isel(pressure_level=level_index)
        result = build_ic.to_model_grid(field).astype(np.float64)
        self._cache[key] = result
        return result

    def surface(self, name: str, when: np.datetime64):
        """Regridded surface truth field, with derived wind10/mslp.

        中文说明：返回插值到模式网格的地表真值场，并派生 wind10/mslp。

        ``name`` is one of wind10, mslp, t2m, ps, tcc, tp.
        中文说明：``name`` 取 wind10、mslp、t2m、ps、tcc、tp 之一。
        """
        key = ("sfc", name, str(when))
        if key in self._cache:
            return self._cache[key]
        dataset, index = self._surface_dataset(when)
        if index < 0:
            return None
        selected = dataset.isel(valid_time=index)
        if name == "wind10":
            # 10 m 风速由 u10、v10 分量合成 / 10 m wind speed from u10 and v10
            u = build_ic.to_model_grid(selected["u10"]).astype(np.float64)
            v = build_ic.to_model_grid(selected["v10"]).astype(np.float64)
            result = np.hypot(u, v)
        elif name == "mslp":
            # Derive MSLP from the same sp/t2m/orography ingredients at every
            # lead, including t0.  Using ERA5's direct msl field at six-hourly
            # leads but a reduction at hourly leads creates an avoidable jump.
            # 中文说明：所有时效（包括 t0）都用相同的 sp/t2m/地形配方推导海平面
            # 气压。若 6 小时时效直接用 ERA5 的 msl、而逐小时时效用静力订正，
            # 会在衔接处产生本可避免的跳变。
            ps = build_ic.to_model_grid(selected["sp"]).astype(np.float64)
            t2m = build_ic.to_model_grid(selected["t2m"]).astype(np.float64)
            result = mslp_reduction(ps, t2m, self.orography())
        else:
            key_name = {"t2m": "t2m", "mslp": "msl", "ps": "sp", "tcc": "tcc", "tp": "tp"}[name]
            result = build_ic.to_model_grid(selected[key_name]).astype(np.float64)
        self._cache[key] = result
        return result

    def orography(self) -> np.ndarray:
        """ERA5 model orography (m); static, taken from the IC surface file.

        中文说明：ERA5 模式地形高度（m）；为静态场，取自初始地表文件。
        """
        z = build_ic.to_model_grid(self.ic_sfc.isel(valid_time=0)["z"]).astype(np.float64)
        return z / G


def mslp_reduction(ps: np.ndarray, t2m: np.ndarray, orography: np.ndarray) -> np.ndarray:
    """Standard synoptic MSLP reduction (both model and ERA5 use their own
    orography and 2 m temperature).

    中文说明：标准的天气学海平面气压订正（模式和 ERA5 各自使用自己的地形与
    2 m 温度）。
    """
    # 先按 6.5 K/km 的温度直减率把 2 m 温度外推到海平面附近（取半地形高度处）。
    # Extrapolate 2 m temperature toward sea level with the 6.5 K/km lapse rate.
    mean_temperature = np.maximum(t2m, 200.0) + 0.0065 * orography / 2.0
    return ps * np.exp(G * orography / (RD * mean_temperature))


# ---------------------------------------------------------------------------
# Model output
# 模式输出
# ---------------------------------------------------------------------------
class Forecast:
    """Frame-stream access layer for one case and initialization method.

    中文说明：某个例、某种初值方法的帧流访问层。
    """

    # 中文说明：检验变量名到帧流变量名的映射。
    # English: mapping from verification names to frame-stream names.
    FRAME_NAMES = {
        "t": "air_temperature",
        "u": "zonal_wind",
        "v": "meridional_wind",
        "q": "specific_humidity",
        "hgt": "geopotential_height",
    }

    def __init__(self, case: common.ForecastCase, method: str):
        self.case = case
        self.method = method
        self.grouped: dict[str, list] = {}
        # 依次读取两段帧流并按时间步排序 / read both segments and sort by step
        for path in [
            common.FRAME_ROOT / case.key / f"{method}_seg1.frames",
            common.FRAME_ROOT / case.key / f"{method}_seg2.frames",
        ]:
            if path.exists():
                for name, frames in frame_io.load_stream(path).items():
                    self.grouped.setdefault(name, []).extend(frames)
        for frames in self.grouped.values():
            frames.sort(key=lambda frame: frame.step)
        self.ic_restart = (
            common.IC_ROOT / case.key / ("direct.restart" if method == "direct" else "nudged.restart")
        )

    def times(self) -> "list[np.datetime64]":
        """Valid times of all archived frames (from air_temperature).

        中文说明：返回所有归档帧的有效时刻（以 air_temperature 帧为准）。
        """
        return [np.datetime64(frame.valid_time) for frame in self.grouped["air_temperature"]]

    def field_at(self, name: str, index: int) -> np.ndarray:
        """Frame field by verification name (or raw stream name).

        中文说明：按检验变量名（或原始帧流名）取第 index 帧的场。
        """
        return self.grouped[self.FRAME_NAMES.get(name, name)][index].values

    def pressure(self, index: int) -> np.ndarray:
        """3-D pressure on the sigma levels of frame ``index``.

        中文说明：第 index 帧对应 σ 层上的三维气压场。
        """
        ps = self.grouped["surface_pressure"][index].values[0]
        return common.SIGMA_FULL[:, None, None] * ps[None, :, :]

    def surface_at(self, name: str, index: int) -> np.ndarray:
        """Surface diagnostic by verification name.

        中文说明：按检验变量名取第 index 帧的地表诊断量。
        """
        if name == "wind10":
            # 10 m 风速用最低模式层风速近似 / approximate 10 m wind with the lowest model level
            u = self.grouped["zonal_wind"][index].values[-1]
            v = self.grouped["meridional_wind"][index].values[-1]
            return np.hypot(u, v).astype(np.float64)
        key = {
            "t2m": "surface_air_temperature",
            "ps": "surface_pressure",
            "tcc": "total_cloud_cover",
            "tp": "precipitation",
        }[name]
        return self.grouped[key][index].values[0].astype(np.float64)


def interpolate_to_pressure(data: np.ndarray, p_model: np.ndarray, level_hpa: float):
    """Linear-in-log-p interpolation of a sigma-level field to one pressure level.

    中文说明：把 σ 层上的场按“对数气压线性”插值到单个气压层。

    Points where the target pressure is below the model surface are masked.
    中文说明：目标气压低于模式地表的点会被掩膜（置 NaN）。
    """
    nlev = data.shape[0]
    p_target = np.full(data.shape[1:], level_hpa * 100.0)
    # 自下而上找到第一层高于目标气压的模式层 / find the first model level above the target
    above = p_model >= p_target[None, :, :]
    idx_lower = np.argmax(above, axis=0)
    idx_upper = np.clip(idx_lower - 1, 0, nlev - 1)
    valid = above.any(axis=0)
    p_lower = np.take_along_axis(p_model, idx_lower[None], axis=0)[0]
    p_upper = np.take_along_axis(p_model, idx_upper[None], axis=0)[0]
    ratio = p_lower / p_upper
    denominator = np.where(ratio > 1.0, np.log(ratio), 1.0)
    # 对数气压线性权重 / linear-in-log-p weight
    w = (np.log(p_target) - np.log(p_upper)) / denominator
    w = np.clip(w, 0.0, 1.0)
    f_lower = np.take_along_axis(data, idx_lower[None], axis=0)[0]
    f_upper = np.take_along_axis(data, idx_upper[None], axis=0)[0]
    result = f_upper + w * (f_lower - f_upper)
    result[~valid] = np.nan
    return result


def hydrostatic_height(t: np.ndarray, q: np.ndarray, ps: np.ndarray,
                       orography: np.ndarray) -> np.ndarray:
    """Geopotential height (m) on the sigma levels of a restart.

    中文说明：按重启文件 σ 层计算位势高度（m）。

    The restart carries no height record, so the initial 500 hPa height map
    is rebuilt with the hypsometric equation from the virtual temperature.
    The lowest layer starts at the model orography; each level above it is
    obtained by integrating ``R*Tv*dln(p)`` between layer centres.

    中文说明：重启文件不含高度记录，因此用虚温通过压高公式重建初始 500 hPa
    高度场。最低层从模式地形起算，其上每一层由层中心之间积分 ``R*Tv*dln(p)``
    得到。
    """
    virtual_temperature = t * (1.0 + 0.608 * q)
    pressure = common.SIGMA_FULL[:, None, None] * ps[None, :, :]
    scale = RD / common.GRAVITY
    height = np.empty_like(t)
    # 最低层：从地形高度积分到该层中心 / lowest layer: integrate from the surface to the layer centre
    height[-1] = orography + scale * virtual_temperature[-1] * np.log(ps / pressure[-1])
    for k in range(common.NLEV - 2, -1, -1):
        mean_temperature = 0.5 * (virtual_temperature[k] + virtual_temperature[k + 1])
        height[k] = height[k + 1] + scale * mean_temperature * np.log(
            pressure[k + 1] / pressure[k]
        )
    return height


def metrics(model: np.ndarray, truth: np.ndarray) -> tuple[float, float]:
    """Area-weighted RMSE and bias over the valid points.

    中文说明：在有效格点上计算面积加权的 RMSE 与偏差。

    Returns ``(nan, nan)`` when either field is missing or no valid point
    remains after masking.
    中文说明：若任一输入缺失、或掩膜后没有有效格点，则返回 ``(nan, nan)``。
    """
    if model is None or truth is None:
        return np.nan, np.nan
    mask = np.isfinite(model) & np.isfinite(truth)
    if mask.sum() == 0:
        return np.nan, np.nan
    difference = model - truth
    # Gaussian rows represent unequal latitude-band areas. Renormalise after
    # masking so missing below-ground values do not alter the denominator.
    # 中文说明：高斯网格各行代表面积不等的纬度带；掩膜后重新归一化，使地下
    # 缺测值不会影响分母。
    weights = build_ic.MODEL_WEIGHTS[:, None]
    valid_weights = np.broadcast_to(weights, difference.shape)
    valid_weights = np.where(mask, valid_weights, 0.0)
    denominator = valid_weights.sum()
    if denominator <= 0.0:
        return np.nan, np.nan
    mse = np.sum(valid_weights * np.where(mask, difference**2, 0.0)) / denominator
    bias = np.sum(valid_weights * np.where(mask, difference, 0.0)) / denominator
    return float(np.sqrt(mse)), float(bias)


# ---------------------------------------------------------------------------
# Initial condition synthesis (lead time 0)
# 初始条件重构（预报时效 0）
# ---------------------------------------------------------------------------
def synth_restart(path: Path) -> dict[str, np.ndarray]:
    """Reconstruct gridded fields from a restart's spectral coefficients.

    中文说明：由重启文件的谱系数重构格点场。

    The restart stores T, divergence, vorticity, humidity and surface pressure
    as spectra; this helper runs the model's synthesis tool and converts the
    output back to physical units (wind speeds are re-dimensionalised with
    ``cv/cos(lat)`` and surface pressure with ``psurf``).
    中文说明：重启文件以谱系数形式保存 T、散度、涡度、比湿和地面气压；本函数
    调用模式的谱合成工具，并把输出换算回物理量（风速乘回 ``cv/cos(纬度)``，
    地面气压乘回 ``psurf``）。
    """
    records = dict(restart_io.read_raw_records(path))
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for name in ("st", "sd", "sz", "sq"):
            values = restart_io.decode_float(records[name], common.NRSP * common.NLEV)
            values.reshape(common.NLEV, common.NRSP).astype("<f4").tofile(tmpdir / f"{name}.bin")
        restart_io.decode_float(records["sp"], common.NRSP).astype("<f4").tofile(tmpdir / "sp.bin")
        subprocess.run(
            [str(common.TOOLS_ROOT / "plasic_era5_tools"), "synth", str(tmpdir)], check=True,
            stdout=subprocess.DEVNULL,
        )
        out = {}
        for key, filename in [("t", "synth_t.bin"), ("u", "synth_u.bin"), ("v", "synth_v.bin"),
                              ("q", "synth_q.bin"), ("ps", "synth_ps.bin")]:
            out[key] = np.fromfile(tmpdir / filename, dtype="<f4").astype(np.float64)
    out["t"] = out["t"].reshape(common.NLEV, common.NLAT, common.NLON)
    out["q"] = out["q"].reshape(common.NLEV, common.NLAT, common.NLON)
    # 无量纲 Robert 形风速 → 物理风速 / dimensionless Robert-form winds -> physical winds
    cos_lat = build_ic.MODEL_COS_LAT[None, :, None]
    out["u"] = out["u"].reshape(common.NLEV, common.NLAT, common.NLON) * common.CV / cos_lat
    out["v"] = out["v"].reshape(common.NLEV, common.NLAT, common.NLON) * common.CV / cos_lat
    # 无量纲地面气压 → Pa / dimensionless surface pressure -> Pa
    psurf = build_ic.model_orography()[1]
    out["ps"] = out["ps"].reshape(common.NLAT, common.NLON) * psurf
    return out


# ---------------------------------------------------------------------------
# Case verification
# 个例检验
# ---------------------------------------------------------------------------
def verify_case(case: common.ForecastCase, method: str, require_exact: bool = False) -> dict:
    """Compute every metric and map for one case/method and return the archive.

    中文说明：计算某个例/初值方法的全部指标与场图，返回归档字典。
    """
    truth = Truth(case, require_exact=require_exact)
    forecast = Forecast(case, method)
    times = forecast.times()
    t0 = np.datetime64(common.init_datetime(case))
    lead_hours = np.array([(t - t0) / np.timedelta64(1, "h") for t in times])
    initial = synth_restart(forecast.ic_restart)
    orography_model = build_ic.model_orography()[0]
    orography_era5 = truth.orography()

    result: dict[str, np.ndarray] = {
        "lead_hours": lead_hours,
        "times": np.array([str(t) for t in times]),
        "metric_area_weighting": np.array("Gaussian quadrature latitude weights; uniform longitude"),
        "verification_temporal_source": np.array(truth.verification_source),
    }
    # 时效 0（初值本身）+ 所有帧时效 / lead 0 (the IC itself) plus all frame leads
    all_leads = np.concatenate([[0.0], lead_hours])
    result["all_lead_hours"] = all_leads

    # --- upper-air fields ---------------------------------------------------
    # --- 高层场 -------------------------------------------------------------
    for name, levels in UPPER_LEVELS.items():
        rmse = np.full((len(levels), all_leads.size), np.nan)
        bias = np.full((len(levels), all_leads.size), np.nan)
        for i, level in enumerate(levels):
            # 时效 0：由重启文件谱合成得到 / lead 0: synthesised from the restart spectra
            model_field = interpolate_to_pressure(
                initial[name], common.SIGMA_FULL[:, None, None] * initial["ps"][None, :, :], level
            )
            rmse[i, 0], bias[i, 0] = metrics(model_field, truth.upper(name, level, t0))
            for j in range(len(times)):
                model_field = interpolate_to_pressure(forecast.field_at(name, j), forecast.pressure(j), level)
                rmse[i, j + 1], bias[i, j + 1] = metrics(model_field, truth.upper(name, level, times[j]))
        result[f"rmse_{name}"] = rmse
        result[f"bias_{name}"] = bias
        result[f"levels_{name}"] = np.array(levels)

    # --- geopotential height at 500 hPa ------------------------------------
    # --- 500 hPa 位势高度 ---------------------------------------------------
    rmse = np.full((1, all_leads.size), np.nan)
    bias = np.full((1, all_leads.size), np.nan)
    rmse[0, 0], bias[0, 0] = np.nan, np.nan  # the IC restart has no height record
                                             # 中文说明：初值重启文件没有高度记录
    for j in range(len(times)):
        model_field = interpolate_to_pressure(forecast.field_at("hgt", j), forecast.pressure(j), 500)
        truth_field = truth.upper("z", 500, times[j])
        if truth_field is not None:
            # ERA5 的 z 是位势，除以 g 得到位势高度（m）。
            # ERA5 z is geopotential; divide by g to obtain geopotential height (m).
            rmse[0, j + 1], bias[0, j + 1] = metrics(model_field, truth_field / G)
    result["rmse_z"] = rmse
    result["bias_z"] = bias
    result["levels_z"] = np.array([500])

    # --- surface fields -----------------------------------------------------
    # --- 地表场 -------------------------------------------------------------
    # 先把所有帧地表场堆叠成 (时刻, 纬度, 经度) 数组，便于统一计算。
    # Stack all surface frames into (time, lat, lon) arrays for uniform handling.
    model_t2m_frames = np.stack([forecast.surface_at("t2m", j) for j in range(len(times))])
    model_ps_frames = np.stack([forecast.surface_at("ps", j) for j in range(len(times))])
    model_wind10_frames = np.stack([forecast.surface_at("wind10", j) for j in range(len(times))])
    model_tcc_frames = np.stack([forecast.surface_at("tcc", j) for j in range(len(times))])
    model_mslp_frames = np.stack(
        [mslp_reduction(model_ps_frames[j], model_t2m_frames[j], orography_model) for j in range(len(times))]
    )
    blank = np.full((common.NLAT, common.NLON), np.nan)
    truth_ps = np.stack([truth.surface("ps", t) if truth.surface("ps", t) is not None else blank
                         for t in times])
    truth_t2m = np.stack([truth.surface("t2m", t) if truth.surface("t2m", t) is not None else blank
                          for t in times])
    truth_wind10 = [truth.surface("wind10", t) for t in times]
    truth_tcc = [truth.surface("tcc", t) for t in times]
    truth_mslp = np.stack(
        [mslp_reduction(truth_ps[j], truth_t2m[j], orography_era5)
         if np.isfinite(truth_ps[j]).any() else blank
         for j in range(len(times))]
    )

    # initial (lead 0) fields: the model's near-surface air temperature is not
    # stored in the restart, so reuse the skin temperature as the best proxy.
    # 中文说明：初值（时效 0）场：重启文件不存近地面气温，因此用皮温作为
    # 最佳近似。
    records = dict(restart_io.read_raw_records(forecast.ic_restart))
    skin = restart_io.decode_float(records["dt"]).reshape(common.NLAT, common.NLON)
    p_initial = common.SIGMA_FULL[:, None, None] * initial["ps"][None, :, :]
    initial_maps = {
        "t500": interpolate_to_pressure(initial["t"], p_initial, 500),
        "t700": interpolate_to_pressure(initial["t"], p_initial, 700),
        "u250": interpolate_to_pressure(initial["u"], p_initial, 250),
        "u850": interpolate_to_pressure(initial["u"], p_initial, 850),
        "v250": interpolate_to_pressure(initial["v"], p_initial, 250),
        "v850": interpolate_to_pressure(initial["v"], p_initial, 850),
        "q700": interpolate_to_pressure(initial["q"], p_initial, 700),
        "mslp": mslp_reduction(initial["ps"], skin, orography_model),
        "t2m": skin,
        # Near-surface wind is represented by the lowest sigma level, the same
        # proxy the frame stream uses for the WeatherBench 2 wind10 metric.
        # 中文说明：近地面风用最低 σ 层表示，与帧流对 WeatherBench 2 wind10
        # 指标所用的近似一致。
        "wind10": np.hypot(initial["u"][-1], initial["v"][-1]),
        # The restart has no height record; rebuild the initial 500 hPa height
        # hydrostatically so the lead-0 map is not blank.
        # 中文说明：重启文件没有高度记录，用静力关系重建初始 500 hPa 高度，
        # 避免时效 0 的场图为空白。
        "z500": interpolate_to_pressure(
            hydrostatic_height(initial["t"], initial["q"], initial["ps"], orography_model),
            p_initial, 500,
        ),
        # The restart has no cloud field; the initial PlaSiC cloud map is blank.
        # 中文说明：重启文件没有云量场，PlaSiC 初始云量图保持空白（NaN）。
        "tcc": np.full((common.NLAT, common.NLON), np.nan),
    }

    surface_pairs = {
        "t2m": (np.concatenate([[skin], model_t2m_frames]),
                np.concatenate([[truth.surface("t2m", t0)], truth_t2m])),
        "mslp": (np.concatenate([[initial_maps["mslp"]], model_mslp_frames]),
                 np.concatenate([[truth.surface("mslp", t0)], truth_mslp])),
        "ps": (np.concatenate([[initial["ps"]], model_ps_frames]),
               np.concatenate([[truth.surface("ps", t0)], truth_ps])),
    }
    for name, (model_values, truth_values) in surface_pairs.items():
        rmse = np.full(all_leads.size, np.nan)
        bias = np.full(all_leads.size, np.nan)
        for j in range(all_leads.size):
            if truth_values[j] is None:
                continue
            rmse[j], bias[j] = metrics(model_values[j], truth_values[j])
        result[f"rmse_{name}"] = rmse
        result[f"bias_{name}"] = bias

    # --- maps at selected lead times ---------------------------------------
    # --- 选定时效的场图 -----------------------------------------------------
    map_frame_indices = [int(np.argmin(np.abs(lead_hours - h))) for h in MAP_TIMES_H]
    # Keep both temperature levels in the archive: t500 remains available to
    # older diagnostics, while the state maps use t700 for a like-for-like
    # comparison with q700.
    # 中文说明：归档同时保留两个温度层：t500 供旧诊断使用，而状态图使用
    # t700，与 q700 形成同层对比。
    for name, level in [("t", 500), ("t", 700), ("u", 250), ("u", 850),
                        ("v", 250), ("v", 850), ("q", 700)]:
        model_maps = [initial_maps[f"{name}{level}"]]
        truth_maps = [truth.upper(name, level, t0)]
        for j in map_frame_indices[1:]:
            model_maps.append(
                interpolate_to_pressure(forecast.field_at(name, j), forecast.pressure(j), level)
            )
            truth_maps.append(truth.upper(name, level, times[j]))
        result[f"map_model_{name}{level}"] = np.stack(model_maps)
        result[f"map_truth_{name}{level}"] = np.stack(
            [np.full((common.NLAT, common.NLON), np.nan) if m is None else m for m in truth_maps]
        )
    result["map_model_z500"] = np.stack(
        [initial_maps["z500"]]
        + [interpolate_to_pressure(forecast.field_at("hgt", j), forecast.pressure(j), 500)
           for j in map_frame_indices[1:]]
    )
    result["map_truth_z500"] = np.stack(
        [
            truth.upper("z", 500, t0) / G if truth.upper("z", 500, t0) is not None
            else np.full((common.NLAT, common.NLON), np.nan)
        ]
        + [
            truth.upper("z", 500, times[j]) / G if truth.upper("z", 500, times[j]) is not None
            else np.full((common.NLAT, common.NLON), np.nan)
            for j in map_frame_indices[1:]
        ]
    )
    result["map_model_t2m"] = np.stack([initial_maps["t2m"]] + [model_t2m_frames[j] for j in map_frame_indices[1:]])
    result["map_truth_t2m"] = np.stack(
        [truth.surface("t2m", t0)]
        + [
            truth_t2m[j] if truth_t2m[j] is not None else np.full((common.NLAT, common.NLON), np.nan)
            for j in map_frame_indices[1:]
        ]
    )
    result["map_model_mslp"] = np.stack(
        [initial_maps["mslp"]] + [model_mslp_frames[j] for j in map_frame_indices[1:]]
    )
    result["map_truth_mslp"] = np.stack(
        [truth.surface("mslp", t0)]
        + [
            truth_mslp[j] if truth_mslp[j] is not None else np.full((common.NLAT, common.NLON), np.nan)
            for j in map_frame_indices[1:]
        ]
    )
    result["map_model_wind10"] = np.stack(
        [initial_maps["wind10"]] + [model_wind10_frames[j] for j in map_frame_indices[1:]]
    )
    result["map_truth_wind10"] = np.stack(
        [truth.surface("wind10", t0) if truth.surface("wind10", t0) is not None
         else np.full((common.NLAT, common.NLON), np.nan)]
        + [
            truth_wind10[j] if truth_wind10[j] is not None else np.full((common.NLAT, common.NLON), np.nan)
            for j in map_frame_indices[1:]
        ]
    )
    result["map_model_tcc"] = np.stack(
        [initial_maps["tcc"]] + [model_tcc_frames[j] for j in map_frame_indices[1:]]
    )
    result["map_truth_tcc"] = np.stack(
        [truth.surface("tcc", t0) if truth.surface("tcc", t0) is not None
         else np.full((common.NLAT, common.NLON), np.nan)]
        + [
            truth_tcc[j] if truth_tcc[j] is not None else np.full((common.NLAT, common.NLON), np.nan)
            for j in map_frame_indices[1:]
        ]
    )
    result["map_times_h"] = np.array(MAP_TIMES_H, dtype=np.float64)

    # --- persistence baseline (initial state kept constant) ------------------
    # --- 持续性基线（初值状态保持不变） --------------------------------------
    persistence_specs = [
        ("t", 850, initial["t"], "t"),
        ("t", 500, initial["t"], "t"),
        ("t", 200, initial["t"], "t"),
        ("u", 250, initial["u"], "u"),
        ("q", 700, initial["q"], "q"),
    ]
    for name, level, source, source_name in persistence_specs:
        model_field = interpolate_to_pressure(source, p_initial, level)
        rmse = np.full(all_leads.size, np.nan)
        for j in range(all_leads.size):
            # 持续性基线的比较对象与目标时刻相差一个时效 / compare persistence with truth one lead earlier
            when = t0 if j == 0 else times[j - 1]
            truth_field = truth.upper(source_name, level, when)
            if truth_field is None:
                continue
            rmse[j] = metrics(model_field, truth_field)[0]
        result[f"persistence_{name}{level}"] = rmse
    rmse = np.full(all_leads.size, np.nan)
    for j in range(all_leads.size):
        when = t0 if j == 0 else times[j - 1]
        truth_field = truth.surface("t2m", when)
        if truth_field is not None:
            rmse[j] = metrics(skin, truth_field)[0]
    result["persistence_t2m"] = rmse
    rmse = np.full(all_leads.size, np.nan)
    for j in range(all_leads.size):
        when = t0 if j == 0 else times[j - 1]
        truth_field = truth.surface("mslp", when)
        if truth_field is not None:
            rmse[j] = metrics(initial_maps["mslp"], truth_field)[0]
    result["persistence_mslp"] = rmse
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--methods", nargs="*", default=["direct", "nudged"])
    parser.add_argument(
        "--require-exact",
        action="store_true",
        help="require complete CDS pressure/surface 6-hour files and hourly shock files",
    )
    args = parser.parse_args()
    VERIF_ROOT.mkdir(parents=True, exist_ok=True)
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        for method in args.methods:
            output = VERIF_ROOT / f"{case.key}_{method}.npz"
            print(f"[{case.key}/{method}] -> {output.name}", flush=True)
            result = verify_case(case, method, require_exact=args.require_exact)
            np.savez_compressed(output, **result)
            print(f"    done, lead times {result['lead_hours'].size + 1}")


if __name__ == "__main__":
    main()
