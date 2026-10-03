#!/usr/bin/env python3
"""Download the ERA5 fields needed by the PlaSiC forecast experiment.

中文说明：下载 PlaSiC 预报试验所需的 ERA5 场。

Three jobs per case:

中文说明：每个个例包含三类下载任务：

* ``nudge``   hourly pressure-level fields over the nudging window
              [t0 - 12 h, t0] (t, u, v, q on 32 levels);
  中文说明：``nudge``——松弛同化窗口 [t0-12 h, t0] 内逐小时的气压层场
  （32 层上的 t、u、v、q）；
* ``verif``   6-hourly pressure-level fields over [t0, t0 + 10 d]
              (t, u, v, q, z on 25 levels);
  中文说明：``verif``——[t0, t0+10 d] 内每 6 小时的气压层场
  （25 层上的 t、u、v、q、z）；
* ``surface`` single-level fields at t0 - 12 h and t0 (initial surface state)
              and 6-hourly over [t0, t0 + 10 d] (verification surface state).
  中文说明：``surface``——t0-12 h 与 t0 的单层场（初始地表状态），以及
  [t0, t0+10 d] 内每 6 小时的场（检验用地表状态）。

Everything is downloaded on a 1.5 deg global grid, which is close to the T85
model resolution and keeps the archive small.  Existing files are skipped.

中文说明：所有数据都下载到 1.5° 的全球网格上，该分辨率接近 T85 模式网格，
同时能让归档体积保持在较小规模。已存在且校验通过的文件会被跳过。
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

DATASET_PL = "reanalysis-era5-pressure-levels"  # 气压层数据集 / pressure-level dataset
DATASET_SL = "reanalysis-era5-single-levels"    # 单层数据集 / single-level dataset

# 中文说明：松弛窗口使用的气压层变量（t、u、v、q）。
# English: pressure-level variables used inside the nudging window (t, u, v, q).
PL_VARIABLES = [
    "temperature",
    "u_component_of_wind",
    "v_component_of_wind",
    "specific_humidity",
]
# 中文说明：检验额外需要位势（z），用于 Z500 检验。
# English: verification additionally needs geopotential (z) for Z500 scores.
PL_VARIABLES_VERIF = PL_VARIABLES + ["geopotential"]

# 中文说明：构建初始地表状态所需的单层变量（土壤温湿、海温、海冰、积雪等）。
# English: single-level variables needed to build the initial surface state
# (soil temperature/moisture, SST, sea ice, snow depth, ...).
SL_VARIABLES_IC = [
    "surface_pressure",
    "geopotential",
    "skin_temperature",
    "sea_surface_temperature",
    "2m_temperature",
    "soil_temperature_level_1",
    "soil_temperature_level_2",
    "soil_temperature_level_3",
    "soil_temperature_level_4",
    "volumetric_soil_water_layer_1",
    "volumetric_soil_water_layer_2",
    "volumetric_soil_water_layer_3",
    "volumetric_soil_water_layer_4",
    "snow_depth",
    "sea_ice_cover",
]
# 中文说明：检验所需的单层变量（含降水、总云量、10 m 风等）。
# English: single-level variables needed for verification (precipitation,
# total cloud cover, 10 m wind, ...).
SL_VARIABLES_VERIF = [
    "2m_temperature",
    "surface_pressure",
    "mean_sea_level_pressure",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "skin_temperature",
    "sea_surface_temperature",
    "sea_ice_cover",
    "total_precipitation",
    "total_cloud_cover",
]

# CDS uses long variable names in requests, while the NetCDF files returned by
# ERA5 normally use their GRIB short names.  Keeping this mapping here lets us
# validate an existing file before deciding that it is safe to skip.
# 中文说明：CDS 请求使用变量的长名称，而 ERA5 返回的 NetCDF 文件通常使用
# GRIB 短名。把映射关系集中在这里，就能在决定“可以跳过”之前校验已有文件。
NETCDF_NAMES = {
    "temperature": "t",
    "u_component_of_wind": "u",
    "v_component_of_wind": "v",
    "specific_humidity": "q",
    "geopotential": "z",
    "surface_pressure": "sp",
    "skin_temperature": "skt",
    "sea_surface_temperature": "sst",
    "2m_temperature": "t2m",
    "soil_temperature_level_1": "stl1",
    "soil_temperature_level_2": "stl2",
    "soil_temperature_level_3": "stl3",
    "soil_temperature_level_4": "stl4",
    "volumetric_soil_water_layer_1": "swvl1",
    "volumetric_soil_water_layer_2": "swvl2",
    "volumetric_soil_water_layer_3": "swvl3",
    "volumetric_soil_water_layer_4": "swvl4",
    "snow_depth": "sd",
    # CDS has used both ``ci`` and ``siconc`` for this field in NetCDF
    # exports; accept either short name when validating a completed file.
    # 中文说明：CDS 在 NetCDF 导出中对该变量既用过 ``ci`` 也用过 ``siconc``，
    # 校验完整文件时两种短名都接受。
    "sea_ice_cover": "ci",
    "mean_sea_level_pressure": "msl",
    "10m_u_component_of_wind": "u10",
    "10m_v_component_of_wind": "v10",
    "total_precipitation": "tp",
    "total_cloud_cover": "tcc",
}

AREA = [90.0, -180.0, -90.0, 180.0]  # CDS 区域顺序 [北, 西, 南, 东] / CDS area order [N, W, S, E]


def request_times(times: "list[dt.datetime]") -> dict:
    """Group datetimes into CDS {year, month, day, time} lists.

    中文说明：把日期时间列表整理为 CDS 请求所需的
    {year, month, day, time} 列表，并去重、排序。
    """
    years, months, days, hours = set(), set(), set(), set()
    for when in times:
        years.add(f"{when.year:04d}")
        months.add(f"{when.month:02d}")
        days.add(f"{when.day:02d}")
        hours.add(f"{when.hour:02d}:00")
    return {
        "year": sorted(years),
        "month": sorted(months),
        "day": sorted(days),
        "time": sorted(hours),
    }


def normalize_download(target: Path) -> None:
    """CDS returns a zip when a request mixes instantaneous and accumulated
    variables; merge the members into a single NetCDF file.

    中文说明：当请求同时包含瞬时量和累积量时，CDS 会返回 zip 压缩包；本函数
    把压缩包内的各个成员合并为单个 NetCDF 文件。
    """
    with target.open("rb") as stream:
        if stream.read(2) != b"PK":
            return  # 不是 zip，无需处理 / not a zip archive, nothing to do
    import tempfile

    import xarray as xr

    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(target) as archive:
            archive.extractall(tmp)
        members = sorted(Path(tmp).rglob("*.nc"))
        if not members:
            raise RuntimeError(f"CDS archive {target} contains no NetCDF member")
        datasets = [xr.open_dataset(member) for member in members]
        try:
            merged = xr.merge(datasets, compat="override")
            # Keep a .nc suffix so xarray/netCDF4 can select its writer.
            # 中文说明：临时文件名保留 .nc 后缀，便于 xarray/netCDF4 选择写出后端。
            temporary = target.with_name(f".{target.name}.normalized-{os.getpid()}.tmp.nc")
            try:
                merged.to_netcdf(temporary)
            finally:
                merged.close()
        finally:
            for dataset in datasets:
                dataset.close()
    os.replace(temporary, target)  # 原子替换原文件 / atomic replacement of the target
    print(f"[norm] {target.name} zip -> NetCDF")


def _time_coordinate(dataset):
    """Name of the time coordinate, supporting CDS ``time``/``valid_time``.

    中文说明：返回时间坐标的名称，兼容 CDS 的 ``time``/``valid_time`` 两种写法。
    """
    for name in ("valid_time", "time"):
        if name in dataset.coords:
            return name
    raise ValueError("NetCDF file has no valid_time or time coordinate")


def _validate_file(path: Path, expected_times=None, required_variables=None,
                   pressure_levels=None) -> None:
    """Validate a completed CDS file before it can be skipped.

    中文说明：在允许跳过已有文件之前，先校验该 CDS 文件确实完整。

    CDS jobs can leave a non-empty partial file after a network failure.  A
    non-empty check is therefore insufficient: verify the time axis, required
    variables, and (for pressure-level requests) the requested levels.

    中文说明：CDS 任务在网络故障后可能留下“非空但不完整”的文件，因此仅检查
    文件非空是不够的：还要校验时间轴、必需变量，以及气压层请求中的层次列表。
    """
    import numpy as np
    import xarray as xr

    with xr.open_dataset(path) as dataset:
        time_name = _time_coordinate(dataset)
        if expected_times is not None:
            actual = np.asarray(dataset[time_name].values).astype("datetime64[ns]")
            expected = np.asarray(
                [np.datetime64(when, "ns") for when in expected_times], dtype="datetime64[ns]"
            )
            if set(actual.tolist()) != set(expected.tolist()):
                raise ValueError(
                    f"time axis mismatch: expected {len(expected)} timestamps, "
                    f"found {len(actual)}"
                )
        if required_variables:
            names = set(dataset.data_vars)
            # 中文说明：海冰变量在 CDS 导出中可能叫 ci 或 siconc，需要同时接受。
            # English: the sea-ice variable may be called ci or siconc.
            aliases = {"sea_ice_cover": ("ci", "siconc")}
            missing = [
                variable for variable in required_variables
                if not any(
                    candidate in names
                    for candidate in aliases.get(
                        variable, (NETCDF_NAMES.get(variable, variable), variable)
                    )
                )
            ]
            if missing:
                raise ValueError(f"missing variables: {', '.join(missing)}")
        if pressure_levels is not None:
            level_name = next(
                (name for name in ("pressure_level", "level") if name in dataset.coords),
                None,
            )
            if level_name is None:
                raise ValueError("pressure-level file has no pressure_level coordinate")
            actual = np.asarray(dataset[level_name].values, dtype=float)
            expected = np.asarray(pressure_levels, dtype=float)
            if set(np.round(actual, 6).tolist()) != set(np.round(expected, 6).tolist()):
                raise ValueError(
                    f"pressure levels mismatch: expected {len(expected)}, found {len(actual)}"
                )


def download(client: cdsapi.Client, dataset: str, request: dict, target: Path,
             force: bool = False, expected_times=None, required_variables=None,
             pressure_levels=None, retries: int = 3) -> None:
    """Retrieve one file atomically and validate its exact contents.

    中文说明：以原子方式获取单个文件，并校验其内容完全符合预期。

    ``cdsapi`` writes directly to the path it is given.  We instead download
    to a hidden sibling and replace the destination only after normalization
    and validation succeed.  This prevents interrupted downloads from being
    mistaken for complete archives on the next run.

    中文说明：``cdsapi`` 原本直接写入给定路径；这里改为先下载到隐藏的临时
    文件，只有归一化和校验都通过后才替换目标文件。这样可以防止中断的下载
    在下次运行时被误认为是完整归档。
    """
    if retries < 1:
        raise ValueError("retries must be at least one")
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0 and not force:
        try:
            # Legacy files may still be zip archives.  Normalize them before
            # validation; corrupt archives fall through to a clean download.
            # 中文说明：旧文件可能仍是 zip 归档，校验前先归一化；损坏的归档会
            # 落到下面重新下载。
            normalize_download(target)
            _validate_file(target, expected_times, required_variables, pressure_levels)
        except Exception as exc:  # noqa: BLE001 - stale files must be replaced
            print(f"[redo ] {target.name}: {exc}", flush=True)
        else:
            print(f"[skip] {target.name} ({target.stat().st_size / 1e6:.1f} MB)")
            return

    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        temporary = target.with_name(
            f".{target.name}.download-{os.getpid()}-{attempt}.tmp"
        )
        try:
            temporary.unlink(missing_ok=True)
            print(f"[get ] {target.name} (attempt {attempt}/{retries}) ...", flush=True)
            client.retrieve(dataset, request, str(temporary))
            if not temporary.exists() or temporary.stat().st_size == 0:
                raise RuntimeError("CDS returned an empty file")
            normalize_download(temporary)
            _validate_file(temporary, expected_times, required_variables, pressure_levels)
            os.replace(temporary, target)  # 校验通过后再原子替换 / atomic replace after validation
            print(f"[done] {target.name} ({target.stat().st_size / 1e6:.1f} MB)")
            return
        except Exception as exc:  # noqa: BLE001 - retry the CDS transaction
            last_error = exc
            temporary.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(2 ** (attempt - 1))  # 指数退避 / exponential backoff
                continue
            raise RuntimeError(f"failed to download {target}: {exc}") from exc
    raise RuntimeError(f"failed to download {target}: {last_error}")


def _daily_times(times: "list[dt.datetime]") -> dict[dt.date, list[dt.datetime]]:
    """Group exact requested timestamps by day.

    中文说明：把请求的精确时刻按“天”分组。

    CDS treats the ``day`` and ``time`` arrays as a Cartesian product.  A
    single request spanning several days therefore silently downloads hours
    that were never requested.  Daily requests are a little more numerous,
    but preserve the exact time axis and make custom cases reproducible.

    中文说明：CDS 把 ``day`` 和 ``time`` 数组当作笛卡尔积处理，因此跨多天的
    单个请求会悄悄下载并未要求的时刻。按天请求虽然请求数略多，但能保持精确
    的时间轴，并使自定义个例可以复现。
    """
    grouped: dict[dt.date, list[dt.datetime]] = {}
    for when in times:
        grouped.setdefault(when.date(), []).append(when)
    return grouped


def _download_exact_days(client: cdsapi.Client, dataset: str, base_request: dict,
                         times: "list[dt.datetime]", target: Path,
                         force: bool = False, required_variables=None,
                         pressure_levels=None) -> None:
    """Download exact timestamps day by day, one CDS request per day.

    中文说明：按天逐个下载精确时刻，每天一个 CDS 请求。

    A single-day request keeps the original file name; multi-day requests are
    split into ``stem_YYYYMMDD.suffix`` files.
    中文说明：单天请求保持原文件名；多天请求拆分为
    ``stem_YYYYMMDD.suffix`` 形式的文件。
    """
    grouped = _daily_times(times)
    if len(grouped) == 1:
        day = next(iter(grouped))
        request = {**base_request, **request_times(grouped[day])}
        download(client, dataset, request, target, force=force,
                 expected_times=grouped[day], required_variables=required_variables,
                 pressure_levels=pressure_levels)
        return
    for day, day_times in sorted(grouped.items()):
        request = {**base_request, **request_times(day_times)}
        download(client, dataset, request,
                 target.with_name(f"{target.stem}_{day:%Y%m%d}{target.suffix}"),
                 force=force, expected_times=day_times,
                 required_variables=required_variables,
                 pressure_levels=pressure_levels)


def run_case(client: cdsapi.Client, case: common.ForecastCase, force: bool = False) -> None:
    """Download every ERA5 product required by one forecast case.

    中文说明：下载单个预报个例所需的全部 ERA5 数据产品。
    """
    t0 = common.init_datetime(case)
    directory = common.ERA5_ROOT / case.key
    directory.mkdir(parents=True, exist_ok=True)

    # --- nudging window: hourly pressure levels -----------------------------
    # --- 松弛同化窗口：逐小时气压层场 ----------------------------------------
    # 窗口包含 t0-12 h 到 t0 共 13 个整点 / the window spans t0-12 h .. t0 (13 hours)
    window = [t0 - dt.timedelta(hours=h) for h in range(common.NUDGE_WINDOW_HOURS, -1, -1)]
    for day in sorted({when.date() for when in window}):
        hours = [when for when in window if when.date() == day]
        request = {
            "product_type": ["reanalysis"],
            "variable": PL_VARIABLES,
            "pressure_level": [str(level) for level in common.ERA5_IC_LEVELS],
            "data_format": "netcdf",
            "download_format": "unarchived",
            "area": AREA,
            "grid": [common.ERA5_GRID, common.ERA5_GRID],
            **request_times(hours),
        }
        download(
            client,
            DATASET_PL,
            request,
            directory / f"nudge_pl_{day:%Y%m%d}.nc",
            expected_times=hours,
            required_variables=PL_VARIABLES,
            pressure_levels=common.ERA5_IC_LEVELS,
            force=force,
        )

    # --- verification: 6-hourly pressure levels ------------------------------
    # --- 检验：每 6 小时的气压层场 -------------------------------------------
    # t0 到 t0+10 d，每 6 小时一个时刻 / 6-hourly times from t0 to t0+10 d
    verif_times = [t0 + dt.timedelta(hours=6 * k) for k in range(0, 4 * common.FORECAST_DAYS + 1)]
    base_request = {
        "product_type": ["reanalysis"],
        "variable": PL_VARIABLES_VERIF,
        "pressure_level": [str(level) for level in common.ERA5_VERIF_LEVELS],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": AREA,
        "grid": [common.ERA5_GRID, common.ERA5_GRID],
    }
    # Keep a distinct stem so an older 12-hour archive cannot be mistaken for
    # the new cadence.  The verifier prefers these files when present.
    # 中文说明：使用独立的文件名，使旧的 12 小时归档不会被误当成新的 6 小时
    # 节奏；检验脚本在存在这些文件时优先使用它们。
    _download_exact_days(
        client,
        DATASET_PL,
        base_request,
        verif_times,
        directory / "verif6_pl.nc",
        required_variables=PL_VARIABLES_VERIF,
        pressure_levels=common.ERA5_VERIF_LEVELS,
        force=force,
    )

    # --- surface initial state ----------------------------------------------
    # --- 地表初始状态 --------------------------------------------------------
    base_request = {
        "product_type": ["reanalysis"],
        "variable": SL_VARIABLES_IC,
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": AREA,
        "grid": [common.ERA5_GRID, common.ERA5_GRID],
    }
    # 只下载窗口起点与 t0 两个时刻 / only the window start and t0 are needed
    _download_exact_days(
        client, DATASET_SL, base_request,
        [t0 - dt.timedelta(hours=common.NUDGE_WINDOW_HOURS), t0],
        directory / "surface_ic.nc",
        required_variables=SL_VARIABLES_IC,
        force=force,
    )

    # --- surface verification -----------------------------------------------
    # --- 地表检验场 ----------------------------------------------------------
    base_request = {
        "product_type": ["reanalysis"],
        "variable": SL_VARIABLES_VERIF,
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": AREA,
        "grid": [common.ERA5_GRID, common.ERA5_GRID],
    }
    _download_exact_days(
        client,
        DATASET_SL,
        base_request,
        verif_times,
        directory / "surface_verif6.nc",
        required_variables=SL_VARIABLES_VERIF,
        force=force,
    )


def main() -> None:
    try:
        import cdsapi
    except ImportError as exc:  # keep local build/verify usable without CDS
        # 中文说明：缺少 cdsapi 时给出友好提示；本地 build/verify 阶段无需 CDS。
        # English: give a friendly message; the local build/verify stages do not need CDS.
        raise SystemExit(
            "The download stage requires the 'cdsapi' package; "
            "build/run/verify stages do not."
        ) from exc
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None,
                        help="case keys (default: all)")
    parser.add_argument("--force", action="store_true",
                        help="redownload valid files as well as stale or incomplete files")
    args = parser.parse_args()

    # Keep transient HTTP errors bounded; the cdsapi default of 500 retries
    # can otherwise sleep for hours after a gateway hiccup.
    # 中文说明：限制瞬时 HTTP 错误的重试次数；cdsapi 默认重试 500 次，网关
    # 抖动后可能会休眠数小时。
    client = cdsapi.Client(retry_max=5, sleep_max=30)
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        print(f"=== case {case.key} ({case.init_time}) ===")
        run_case(client, case, force=args.force)


if __name__ == "__main__":
    main()
