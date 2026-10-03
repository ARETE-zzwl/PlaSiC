#!/usr/bin/env python3
"""WeatherBench 2 headline scores for the PlaSiC ERA5-initialized forecasts.

中文说明：针对 ERA5 初始化的 PlaSiC 预报计算 WeatherBench 2 综合评分。

The eight variables follow the WeatherBench 2 deterministic scorecards
(Rasp et al. 2024, JAMES 16, e2023MS004019): 500 hPa geopotential, 850 hPa
temperature, 700 hPa specific humidity, 850 hPa wind vector, 2 m temperature,
surface pressure, 10 m wind speed and 24 h precipitation (RMSE and SEEPS).

中文说明：八个变量与 WeatherBench 2 确定性预报记分卡一致（Rasp 等，2024，
JAMES 16，e2023MS004019）：500 hPa 位势、850 hPa 温度、700 hPa 比湿、
850 hPa 风矢量、2 m 温度、地面气压、10 m 风速以及 24 小时降水（RMSE 与
SEEPS 两种指标）。

Model fields come from the PlaSiC frame streams, truth from the exact CDS
ERA5 archives already used by ``verify.py``.  The 24 h precipitation truth and
the SEEPS climatology come from the cached WB2 products written by
``fetch_wb2_truth.py``.

中文说明：模式场来自 PlaSiC 帧流，真值来自 ``verify.py`` 已在使用的精确
CDS ERA5 归档；24 小时降水真值与 SEEPS 气候态来自 ``fetch_wb2_truth.py``
写出的 WB2 缓存产品。

Scores are written to ``data/verification/wb2_scores.json`` for the report
tables (``make_report_tables.py``).  Every metric is computed for every
6-hourly lead of every archived case and method.

中文说明：评分结果写入 ``data/verification/wb2_scores.json``，供报告表格
（``make_report_tables.py``）使用。每个已归档个例和初值方法的所有 6 小时
时效都会计算全部指标。

Definitions
-----------
中文说明：定义
-----------
* Area weighting follows the local verifier (Gaussian quadrature latitudes).
  中文说明：面积加权与本地检验脚本一致（纬度方向高斯求积权重）。
* Z500 is the geopotential (g * geopotential height) in m2 s-2, matching WB2.
  中文说明：Z500 为位势（g × 位势高度），单位 m2 s-2，与 WB2 一致。
* The 850 hPa wind vector score is sqrt(area_mean((u_f-u_t)^2+(v_f-v_t)^2)).
  中文说明：850 hPa 风矢量评分为 sqrt(面积平均((u_f-u_t)²+(v_f-v_t)²))。
* 10 m wind speed uses the lowest model level as the surface-wind proxy;
  PlaSiC has no dedicated 10 m diagnostic.
  中文说明：10 m 风速用最低模式层作为近地面风近似；PlaSiC 没有专门的
  10 m 诊断量。
* The 24 h model precipitation is the trapezoidal time integral of the
  precipitation-rate frames over (valid - 24 h, valid]; the diagnostic is an
  instantaneous mm/hr rate, so the 6-hourly sampling of days 2-10 makes this
  an approximation.
  中文说明：模式的 24 小时降水由降水率帧在 (valid-24 h, valid] 区间上做梯形
  时间积分得到；该诊断量是瞬时 mm/hr 速率，因此第 2–10 天的 6 小时采样使其
  只是一个近似。
* SEEPS follows the WeatherBench 2 implementation: dry threshold 0.25 mm,
  climatological wet threshold and dry fraction from the 1990-2019 ERA5
  climatology, points with dry fraction outside [0.1, 0.85] masked out.
  中文说明：SEEPS 遵循 WeatherBench 2 的实现：干阈值 0.25 mm，湿阈值与干
  比例取自 1990–2019 年 ERA5 气候态，干比例在 [0.1, 0.85] 之外的格点被
  掩膜剔除。
"""
from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_ic  # noqa: E402
import common  # noqa: E402
import frame_io  # noqa: E402
import verify  # noqa: E402

TRUTH_ROOT = common.DATA_ROOT / "wb2_truth"                 # WB2 真值缓存目录 / cached WB2 truth
OUTPUT = common.DATA_ROOT / "verification" / "wb2_scores.json"  # 评分输出 / score output
LEADS = np.arange(6, 240 + 6, 6)                             # 评分时效（6 小时间隔）/ scored leads
# 中文说明：评分只需要这些帧流变量，用于跳过其余数据体以节省内存。
# English: only these frame variables are needed for scoring; other payloads are skipped.
WANTED_FRAMES = {
    "air_temperature",
    "specific_humidity",
    "zonal_wind",
    "meridional_wind",
    "geopotential_height",
    "surface_pressure",
    "surface_air_temperature",
    "precipitation",
}
VARIABLES = ["z500", "t850", "q700", "wind850", "t2m", "sp", "wind10",
             "precip24_rmse", "precip24_seeps"]
UNITS = {
    "z500": "m2 s-2",
    "t850": "K",
    "q700": "g kg-1",
    "wind850": "m s-1",
    "t2m": "K",
    "sp": "Pa",
    "wind10": "m s-1",
    "precip24_rmse": "mm",
    "precip24_seeps": "dimensionless",
}
G = common.GRAVITY  # 重力加速度 / gravitational acceleration


# ---------------------------------------------------------------------------
# Frame reading
# 帧流读取
# ---------------------------------------------------------------------------
def iter_selected_frames(path: Path, wanted: set[str]):
    """Yield frames of the wanted variables, skipping the other payloads.

    中文说明：只生成所需变量的帧，其余数据体通过 seek 跳过以节省时间与内存。
    """
    with Path(path).open("rb") as stream:
        header = stream.read(frame_io.FILE_HEADER_BYTES)
        if len(header) != frame_io.FILE_HEADER_BYTES or header[:8] != b"PLASICPS":
            raise ValueError(f"{path} is not a PlaSiC frame stream")
        while True:
            frame_header = stream.read(frame_io.FRAME_HEADER_BYTES)
            if frame_header == b"":
                break  # 正常读到文件尾 / clean end of stream
            if (len(frame_header) != frame_io.FRAME_HEADER_BYTES
                    or frame_header[:8] != b"PFRAME1\0"):
                raise ValueError("corrupt frame header")
            # 注意：这里只解出 payload_bytes（偏移 16），无需读取变量编号。
            # Note: only payload_bytes is unpacked (offset 16); the variable id is not needed.
            (payload_bytes,) = struct.unpack_from("<I", frame_header, 16)
            (step,) = struct.unpack_from("<Q", frame_header, 20)
            date = struct.unpack_from("<IIIII", frame_header, 28)
            (_, frame_nlon, frame_nlat, layers) = struct.unpack_from(
                "<IIII", frame_header, 48)
            name = frame_header[76:108].split(b"\0")[0].decode("ascii", "replace")
            unit = frame_header[108:124].split(b"\0")[0].decode("ascii", "replace")
            if name not in wanted:
                stream.seek(payload_bytes, 1)  # 跳过无关变量 / skip unwanted payload
                continue
            payload = stream.read(payload_bytes)
            if len(payload) != payload_bytes:
                raise ValueError("truncated frame payload")
            values = np.frombuffer(payload, dtype="<f4").reshape(
                layers, frame_nlat, frame_nlon).copy()
            yield frame_io.Frame(
                variable_id=0, name=name, unit=unit, step=step,
                date=tuple(int(v) for v in date), level_index=0, layers=layers,
                nlon=frame_nlon, nlat=frame_nlat, values=values,
            )


def load_frames(case: common.ForecastCase, method: str) -> dict[str, list]:
    """Load both frame segments and sort each variable by model step.

    中文说明：读取两段帧流，并按模式时间步对每个变量排序。
    """
    grouped: dict[str, list] = {}
    for path in [
        common.FRAME_ROOT / case.key / f"{method}_seg1.frames",
        common.FRAME_ROOT / case.key / f"{method}_seg2.frames",
    ]:
        if not path.exists():
            continue
        for frame in iter_selected_frames(path, WANTED_FRAMES):
            grouped.setdefault(frame.name, []).append(frame)
    for frames in grouped.values():
        frames.sort(key=lambda frame: frame.step)
    return grouped


def frame_times_hours(case: common.ForecastCase, frames: dict) -> np.ndarray:
    """Lead hours of every loaded frame relative to t0.

    中文说明：返回每个已读取帧相对 t0 的预报时效（小时）。
    """
    times = [np.datetime64(frame.valid_time) for frame in frames["air_temperature"]]
    t0 = np.datetime64(common.init_datetime(case))
    return np.asarray([(t - t0) / np.timedelta64(1, "h") for t in times],
                      dtype=np.float64)


# ---------------------------------------------------------------------------
# Fields
# 场构造
# ---------------------------------------------------------------------------
def model_fields(frames: dict, index: int) -> dict[str, np.ndarray]:
    """All WB2 model diagnostics at frame ``index``.

    中文说明：第 index 帧对应的全部 WB2 模式诊断量。
    """
    ps = frames["surface_pressure"][index].values[0].astype(np.float64)
    # 由 σ 层与地面气压得到三维气压 / 3-D pressure from sigma levels and surface pressure
    pressure = common.SIGMA_FULL[:, None, None] * ps[None, :, :]
    return {
        # z500 保持位势（× g），与 WB2 定义一致 / keep z500 as geopotential (x g)
        "z500": G * verify.interpolate_to_pressure(
            frames["geopotential_height"][index].values, pressure, 500),
        "t850": verify.interpolate_to_pressure(
            frames["air_temperature"][index].values, pressure, 850),
        # 比湿 kg/kg → g/kg / specific humidity kg/kg -> g/kg
        "q700": verify.interpolate_to_pressure(
            frames["specific_humidity"][index].values, pressure, 700) * 1000.0,
        "u850": verify.interpolate_to_pressure(
            frames["zonal_wind"][index].values, pressure, 850),
        "v850": verify.interpolate_to_pressure(
            frames["meridional_wind"][index].values, pressure, 850),
        "t2m": frames["surface_air_temperature"][index].values[0].astype(np.float64),
        "sp": ps,
        # 10 m 风速用最低模式层风速近似 / approximate 10 m wind with the lowest model level
        "wind10": np.hypot(
            frames["zonal_wind"][index].values[-1],
            frames["meridional_wind"][index].values[-1]).astype(np.float64),
    }


def truth_fields(truth: verify.Truth, when: np.datetime64) -> dict[str, np.ndarray]:
    """All WB2 truth diagnostics at valid time ``when``.

    中文说明：``when`` 时刻对应的全部 WB2 真值诊断量。
    """
    return {
        "z500": truth.upper("z", 500, when),
        "t850": truth.upper("t", 850, when),
        "q700": truth.upper("q", 700, when) * 1000.0,  # kg/kg → g/kg
        "u850": truth.upper("u", 850, when),
        "v850": truth.upper("v", 850, when),
        "t2m": truth.surface("t2m", when),
        "sp": truth.surface("ps", when),
        "wind10": truth.surface("wind10", when),
    }


def area_metrics(model, truth):
    """Gaussian-weighted RMSE and bias, identical to verify.metrics.

    中文说明：高斯加权的 RMSE 与偏差，与 verify.metrics 完全一致。
    """
    return verify.metrics(model, truth)


def wind_vector_rmse(model_u, model_v, truth_u, truth_v) -> float:
    """Area-weighted RMSE of the horizontal wind vector.

    中文说明：水平风矢量的面积加权 RMSE。

    Unlike two scalar RMSEs, this penalises the vector difference
    ``(du^2 + dv^2)`` directly, as required by the WB2 wind-vector score.
    中文说明：与分别计算两个标量 RMSE 不同，这里直接惩罚矢量差
    ``(du² + dv²)``，符合 WB2 风矢量评分的要求。
    """
    mask = (np.isfinite(model_u) & np.isfinite(model_v)
            & np.isfinite(truth_u) & np.isfinite(truth_v))
    if not mask.any():
        return float("nan")
    weights = np.broadcast_to(build_ic.MODEL_WEIGHTS[:, None], model_u.shape)
    weights = np.where(mask, weights, 0.0)
    denominator = weights.sum()
    squared = np.where(mask, (model_u - truth_u) ** 2 + (model_v - truth_v) ** 2, 0.0)
    return float(np.sqrt(np.sum(weights * squared) / denominator))


# ---------------------------------------------------------------------------
# 24-hour precipitation and SEEPS
# 24 小时降水与 SEEPS
# ---------------------------------------------------------------------------
def precipitation_integrals(frames: dict, hours: np.ndarray) -> np.ndarray:
    """Cumulative trapezoidal integral of the precipitation rate (mm).

    中文说明：对降水率做梯形法累积积分，得到累积降水量（mm）。
    """
    rates = np.stack([frame.values[0] for frame in frames["precipitation"]])
    rates = rates.astype(np.float64)
    deltas = np.diff(hours)[:, None, None]
    # 相邻帧降水率均值 × 时间间隔 = 区间增量 / mean rate x interval = increment
    increments = 0.5 * (rates[1:] + rates[:-1]) * deltas
    return np.concatenate([np.zeros_like(rates[:1]), np.cumsum(increments, axis=0)])


def window_accumulation(integrals: np.ndarray, hours: np.ndarray,
                        index: int, window_hours: float) -> np.ndarray:
    """24 h accumulation ending at frame ``index`` (left edge extended).

    中文说明：以第 index 帧为终点、向前 24 小时的累积降水（左端点做了延伸）。
    """
    end = hours[index]
    start = end - window_hours
    if start <= hours[0]:
        # Left rectangle using the first available rate for the missing part.
        # 中文说明：缺失的起始时段用第一个可用降水率做左矩形近似。
        missing = hours[0] - start
        # Estimate the pre-first-frame depth from the first increment slope.
        # 中文说明：用第一个区间的增量斜率估计第一帧之前的降水量。
        slope = (integrals[1] - integrals[0]) / (hours[1] - hours[0])
        lead_in = slope * missing
        return integrals[index] + lead_in
    # 在积分曲线上插出左端点值，再相减得到窗口增量。
    # Interpolate the left edge on the integral curve and subtract.
    upper = int(np.searchsorted(hours, start, side="right"))
    lower = upper - 1
    weight = (start - hours[lower]) / (hours[upper] - hours[lower])
    start_value = (1.0 - weight) * integrals[lower] + weight * integrals[upper]
    return integrals[index] - start_value


def seeps_score(model_mm: np.ndarray, truth_mm: np.ndarray,
                threshold_mm: np.ndarray, p1: np.ndarray,
                dry_mm: float = 0.25) -> float:
    """Spatially averaged SEEPS (WeatherBench 2 implementation).

    中文说明：空间平均的 SEEPS 评分（WeatherBench 2 实现）。

    Both forecast and truth are classified into dry / light / heavy categories
    and scored against the climatological dry fraction ``p1``.  Only points
    with 0.1 < p1 < 0.85 contribute.
    中文说明：预报与真值都被划分为干/小雨/大雨三类，并对照气候干比例 ``p1``
    评分；只有 0.1 < p1 < 0.85 的格点参与平均。
    """
    valid = (np.isfinite(model_mm) & np.isfinite(truth_mm)
             & (p1 > 0.1) & (p1 < 0.85))
    if not valid.any():
        return float("nan")
    # 三类掩膜：干（<0.25 mm）、小雨（阈值以下）、大雨（阈值以上）。
    # Three category masks: dry, light (below threshold) and heavy (above threshold).
    f_dry = model_mm < dry_mm
    f_light = (model_mm > dry_mm) & (model_mm < threshold_mm)
    f_heavy = model_mm >= threshold_mm
    t_dry = truth_mm < dry_mm
    t_light = (truth_mm > dry_mm) & (truth_mm < threshold_mm)
    t_heavy = truth_mm >= threshold_mm
    # WB2 的 SEEPS 罚分表；用 errstate 抑制 p1=0/1 边界处的除零警告。
    # WB2 SEEPS penalty table; errstate suppresses divide-by-zero at p1 boundaries.
    with np.errstate(divide="ignore", invalid="ignore"):
        score = 0.5 * (
            f_dry * (t_light * (1.0 / (1.0 - p1)) + t_heavy * (4.0 / (1.0 - p1)))
            + f_light * (t_dry * (1.0 / p1) + t_heavy * (3.0 / (1.0 - p1)))
            + f_heavy * (t_dry * (1.0 / p1 + 3.0 / (2.0 + p1))
                         + t_light * (3.0 / (2.0 + p1)))
        )
    # 纬度高斯加权后做空间平均 / spatial average with Gaussian latitude weights
    weights = np.where(valid, build_ic.MODEL_WEIGHTS[:, None], 0.0)
    return float(np.sum(weights * np.where(valid, score, 0.0)) / weights.sum())


# ---------------------------------------------------------------------------
# Per case / method evaluation
# 逐个例/方法评估
# ---------------------------------------------------------------------------
def measure(case: common.ForecastCase, method: str, truth: verify.Truth,
            cached: dict) -> dict[str, list[float]]:
    """Compute every WB2 score at every 6-hourly lead for one case/method.

    中文说明：计算某个例/初值方法在每个 6 小时时效上的全部 WB2 评分。
    """
    frames = load_frames(case, method)
    if not frames.get("air_temperature"):
        raise RuntimeError(f"no frames for {case.key}/{method}")
    hours = frame_times_hours(case, frames)
    t0 = np.datetime64(common.init_datetime(case))
    integrals = precipitation_integrals(frames, hours)
    precip_truth = cached["precip24_truth_mm"]
    thresholds = cached["precip24_threshold_mm"]
    p1 = cached["p1"]

    results = {name: [] for name in VARIABLES}
    for index, hour in enumerate(LEADS):
        # 找到最接近目标时效的帧；偏差超过 0.5 h 视为缺测。
        # Pick the frame closest to the target lead; > 0.5 h off means missing.
        frame_index = int(np.argmin(np.abs(hours - hour)))
        if abs(hours[frame_index] - hour) > 0.5:
            for name in VARIABLES:
                results[name].append(float("nan"))
            continue
        when = t0 + np.timedelta64(int(hour), "h")
        model = model_fields(frames, frame_index)
        reference = truth_fields(truth, when)
        rmse = {
            name: area_metrics(model[name], reference[name])[0]
            for name in ("z500", "t850", "q700", "t2m", "sp", "wind10")
        }
        rmse["wind850"] = wind_vector_rmse(
            model["u850"], model["v850"], reference["u850"], reference["v850"])
        model_precip = window_accumulation(integrals, hours, frame_index, 24.0)
        truth_precip = precip_truth[index]
        # 降水平方误差同样按纬度高斯加权 / precipitation RMSE uses the same area weights
        rmse["precip24_rmse"] = float(np.sqrt(np.average(
            (model_precip - truth_precip) ** 2,
            weights=np.broadcast_to(build_ic.MODEL_WEIGHTS[:, None],
                                    model_precip.shape))))
        rmse["precip24_seeps"] = seeps_score(
            model_precip, truth_precip, thresholds[index], p1)
        for name in VARIABLES:
            results[name].append(rmse[name])
        print(f"    {case.key}/{method} +{int(hour)} h: "
              f"z500 {rmse['z500']:.1f}  t850 {rmse['t850']:.2f}  "
              f"seeps {rmse['precip24_seeps']:.3f}", flush=True)
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--methods", nargs="*", default=["direct", "nudged"])
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    # p1 与个例无关，读取一次即可 / p1 is case-independent and loaded once
    p1_cache = np.load(TRUTH_ROOT / "seeps_p1.npz")["p1"]
    report = {
        "leads_hours": [int(h) for h in LEADS],
        "units": UNITS,
        "notes": {
            "wind10": "lowest model level used as the surface-wind proxy",
            "precip24_rmse": "24 h accumulation built by trapezoidal "
                             "integration of the mm/hr rate frames",
        },
        "cases": {},
    }
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        # 严格模式：WB2 评分要求精确的 CDS 6 小时/逐小时归档。
        # Strict mode: WB2 scoring requires the exact CDS 6-hourly/hourly archives.
        truth = verify.Truth(case, require_exact=True)
        cached_file = np.load(TRUTH_ROOT / f"{case.key}_precip24.npz")
        cached = {
            "precip24_truth_mm": cached_file["precip24_truth_mm"],
            "precip24_threshold_mm": cached_file["precip24_threshold_mm"],
            "p1": p1_cache,
        }
        report["cases"][case.key] = {}
        for method in args.methods:
            print(f"[score] {case.key}/{method}", flush=True)
            report["cases"][case.key][method] = measure(case, method, truth, cached)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1))
    print(f"[score] wrote {args.output}")


if __name__ == "__main__":
    main()
