#!/usr/bin/env python3
"""Cache the WeatherBench 2 ERA5 truth products used by ``score_wb2.py``.

中文说明：缓存 ``score_wb2.py`` 所需的 WeatherBench 2（WB2）ERA5 真值产品。

The verifier compares PlaSiC with the same ERA5 files it always uses.  Two
extras are needed for the WB2 headline scores that the local archive does not
carry:

中文说明：检验脚本仍使用与以往相同的 ERA5 文件进行对比；但 WB2 综合评分还
需要本地归档中没有的两类数据：

* the 24-hour precipitation accumulation at each valid time, and
* the SEEPS climatological thresholds (1990-2019 ERA5).

中文说明：
* 每个有效时刻的 24 小时累积降水；
* SEEPS 气候阈值（1990–2019 年 ERA5 气候态）。

Both are published in Google's public WeatherBench 2 bucket.  This script
downloads just the required times, interpolates everything onto the PlaSiC
Gaussian grid and stores small ``.npz`` files under ``data/wb2_truth/`` so the
scoring stage runs offline and reproducibly.

中文说明：这两类数据都发布在 Google 的公开 WeatherBench 2 存储桶中。本脚本
只下载所需的时刻，把所有场插值到 PlaSiC 高斯网格，并在 ``data/wb2_truth/``
下保存为小型 ``.npz`` 文件，使评分阶段可以离线、可复现地运行。

Sources
-------
中文说明：数据来源
---------------
* ``gs://weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr``
  (``total_precipitation_24hr``, derived from ERA5 by the WB2 authors).
  中文说明：（``total_precipitation_24hr``，由 WB2 作者从 ERA5 导出。）
* ``gs://weatherbench2/datasets/era5-hourly-climatology/1990-2019_6h_240x121_equiangular_with_poles_conservative.zarr``
  (``total_precipitation_24hr_seeps_threshold``, ``..._dry_fraction``).
  中文说明：（``total_precipitation_24hr_seeps_threshold`` 与 ``..._dry_fraction``。）

See Rasp et al. (2024), *WeatherBench 2*, JAMES 16, e2023MS004019, and the
method description on the WeatherBench 2 website.

中文说明：参见 Rasp 等（2024），*WeatherBench 2*，JAMES 16，e2023MS004019，
以及 WeatherBench 2 网站上的方法说明。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_ic  # noqa: E402
import common  # noqa: E402

OUTPUT_ROOT = common.DATA_ROOT / "wb2_truth"  # 缓存输出目录 / cache output directory
ERA5_6H = (
    "gs://weatherbench2/datasets/era5/"
    "1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
)  # 6 小时分辨率 ERA5 Zarr（WB2 真值）/ 6-hourly ERA5 Zarr (WB2 truth)
CLIMATOLOGY = (
    "gs://weatherbench2/datasets/era5-hourly-climatology/"
    "1990-2019_6h_240x121_equiangular_with_poles_conservative.zarr"
)  # 1990–2019 气候态 Zarr / 1990-2019 climatology Zarr
STORAGE_OPTIONS = {"token": "anon"}  # 匿名访问公开存储桶 / anonymous access to the public bucket
LEADS_HOURS = np.arange(6, 240 + 6, 6)  # 6 小时间隔的预报时效 / 6-hourly lead times


def to_model_grid(field: xr.DataArray) -> np.ndarray:
    """Bilinear map one (latitude, longitude) field to the T85 grid.

    中文说明：把单个 (纬度, 经度) 场双线性插值到 T85 模式网格。
    """
    if field.dims != ("latitude", "longitude"):
        field = field.transpose("latitude", "longitude")
    return build_ic.to_model_grid(field).astype(np.float32)


def valid_times(case: common.ForecastCase) -> np.ndarray:
    """Exact valid timestamps for all archived lead hours of a case.

    中文说明：返回某个例所有归档预报时效对应的精确有效时刻。
    """
    t0 = np.datetime64(common.init_datetime(case), "ns")
    return np.asarray([t0 + np.timedelta64(int(h), "h") for h in LEADS_HOURS])


def fetch_case(case: common.ForecastCase, era5: xr.Dataset,
               climatology: xr.Dataset) -> None:
    """Extract one case's precipitation truth and SEEPS thresholds.

    中文说明：提取一个个例的降水真值和 SEEPS 阈值。
    """
    times = valid_times(case)
    print(f"[fetch] {case.key}: {len(times)} valid times", flush=True)
    # WB2 提供的 24 小时累积降水在时间维上直接按有效时刻选取。
    # WB2's 24 h precipitation accumulation is selected directly by valid time.
    precip = era5["total_precipitation_24hr"].sel(time=times).load()
    threshold = climatology["total_precipitation_24hr_seeps_threshold"]
    # 不同版本的气候态文件用 ``hour`` 或 ``hour_of_day`` 表示一天中的小时。
    # Different climatology versions name the hour axis ``hour`` or ``hour_of_day``.
    hour_name = "hour" if "hour" in threshold.dims else "hour_of_day"
    hours = np.asarray(threshold[hour_name].values, dtype="int64")
    print(f"    climatology hour axis: {hours}", flush=True)
    fields = []
    for index, when in enumerate(times):
        # 按“小时 + 年内日序”索引气候态阈值，再插值到模式网格。
        # Index the climatological threshold by hour + day-of-year, then regrid.
        hour = int((when.astype("datetime64[h]").astype(int)) % 24)
        dayofyear = int(
            (when.astype("datetime64[D]") - when.astype("datetime64[Y]"))
            .astype("timedelta64[D]").astype(int)
        ) + 1
        selected = threshold.sel({hour_name: hour, "dayofyear": dayofyear}).load()
        fields.append(to_model_grid(selected) * 1000.0)  # m -> mm / 米换算为毫米
    out = {
        "leads_hours": LEADS_HOURS,
        "valid_times": np.asarray([str(np.datetime64(t)) for t in times]),
        "precip24_truth_mm": np.stack([to_model_grid(f) for f in precip]) * 1000.0,
        "precip24_threshold_mm": np.stack(fields),
    }
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUTPUT_ROOT / f"{case.key}_precip24.npz", **out)
    print(f"    wrote {case.key}_precip24.npz", flush=True)


def fetch_dry_fraction(climatology: xr.Dataset) -> None:
    """Time-mean SEEPS dry fraction p1, as used by WeatherBench 2.

    中文说明：计算 SEEPS 评分使用的气候平均干（旱）比例 p1，与
    WeatherBench 2 的实现一致。
    """
    dry = climatology["total_precipitation_24hr_seeps_dry_fraction"]
    # 对一天中的小时和年内日序取平均，得到时间平均的 p1。
    # Average over hour-of-day and day-of-year to obtain the time-mean p1.
    mean = dry.mean([d for d in ("hour", "hour_of_day", "dayofyear") if d in dry.dims])
    p1 = to_model_grid(mean.load())
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUTPUT_ROOT / "seeps_p1.npz", p1=p1)
    print(f"[fetch] wrote seeps_p1.npz (p1 range {p1.min():.3f}-{p1.max():.3f})",
          flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None)
    args = parser.parse_args()
    # 直接打开远程 Zarr（只读、匿名），用完显式关闭。
    # Open the remote Zarr stores read-only and anonymously; close them explicitly.
    era5 = xr.open_zarr(ERA5_6H, storage_options=STORAGE_OPTIONS,
                        consolidated=True)
    climatology = xr.open_zarr(CLIMATOLOGY, storage_options=STORAGE_OPTIONS,
                               consolidated=True)
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        fetch_case(case, era5, climatology)
    # p1 与个例无关，只需抓取一次 / p1 is case-independent and fetched once.
    fetch_dry_fraction(climatology)
    era5.close()
    climatology.close()


if __name__ == "__main__":
    main()
