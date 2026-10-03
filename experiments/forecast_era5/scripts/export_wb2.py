#!/usr/bin/env python3
"""Export a PlaSiC case in the WeatherBench 2 by-init layout.

中文说明：按 WeatherBench 2 的 “by-init” 布局导出单个 PlaSiC 个例。

The export is deliberately small and contains the headline deterministic
fields (Z500, T850 and T2M).  It uses the same 1.5-degree target grid as the
public WB2 ERA5 archive.  The interpolation here is bilinear; a strict WB2
submission should use the WB2 conservative-regridding utility before scoring.

中文说明：导出内容刻意保持精简，只包含核心确定性预报场（Z500、T850 和
T2M）。目标网格与公开的 WB2 ERA5 归档一致，均为 1.5°。这里的插值是双线性
的；若要严格提交 WB2 评分，应先用 WB2 的守恒重网格化工具处理。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import verify  # noqa: E402


REG_LAT = np.arange(-90.0, 90.0 + 1.5, 1.5)  # WB2 目标纬度 / WB2 target latitudes
REG_LON = np.arange(0.0, 360.0, 1.5)         # WB2 目标经度 / WB2 target longitudes


def regular_grid(values: np.ndarray) -> np.ndarray:
    """Bilinearly map [lat, lon] Gaussian fields to the WB2 1.5° grid.

    中文说明：把 [纬, 经] 高斯网格场双线性插值到 WB2 的 1.5° 网格。
    """
    da = xr.DataArray(values, dims=("latitude", "longitude"),
                      coords={"latitude": verify.build_ic.MODEL_LAT,
                              "longitude": verify.build_ic.MODEL_LON})
    # xarray needs monotonic coordinates; keep the periodic seam explicit.
    # 中文说明：xarray 要求坐标单调；这里显式补上 360° 的周期接缝。
    da = da.sortby("latitude").sortby("longitude")
    seam = da.isel(longitude=0).assign_coords(longitude=360.0)
    da = xr.concat([da, seam], dim="longitude")
    return da.interp(latitude=REG_LAT, longitude=REG_LON).values


def export_case(case: common.ForecastCase, method: str, output: Path) -> None:
    """Export one case/method as WB2 forecast + observation Zarr/NetCDF.

    中文说明：把某个例/初值方法导出为 WB2 预报与观测两种 Zarr/NetCDF 文件。
    """
    forecast = verify.Forecast(case, method)
    truth = verify.Truth(case)
    t0 = np.datetime64(common.init_datetime(case))
    leads = np.asarray([(np.datetime64(t) - t0) / np.timedelta64(1, "h")
                        for t in forecast.times()], dtype=np.int64)
    # WB2 by-init products use a regular 6-hour cadence.  Keep the richer
    # hourly first-day stream for local diagnostics, but select the 6-hour
    # subset for this export.
    # 中文说明：WB2 的 by-init 产品采用规则 6 小时间隔。本地诊断保留更丰富的
    # 首日逐小时流，但本导出只选取 6 小时的子集。
    keep = (leads % 6 == 0)
    if leads.size == 0:
        raise RuntimeError(f"no frames found for {case.key}/{method}")

    z500, t850, t2m = [], [], []
    obs_z500, obs_t850, obs_t2m, selected_leads = [], [], [], []
    for j, when in enumerate(forecast.times()):
        if not keep[j]:
            continue
        pressure = forecast.pressure(j)
        truth_z = truth.upper("z", 500, when)
        truth_t = truth.upper("t", 850, when)
        truth_t2 = truth.surface("t2m", when)
        # 三个真值场缺一不可，否则跳过该时刻 / skip a time unless all three truth fields exist
        if truth_z is None or truth_t is None or truth_t2 is None:
            continue
        z500.append(regular_grid(verify.interpolate_to_pressure(
            forecast.field_at("hgt", j), pressure, 500) ))
        t850.append(regular_grid(verify.interpolate_to_pressure(
            forecast.field_at("t", j), pressure, 850) ))
        t2m.append(regular_grid(forecast.surface_at("t2m", j)))
        # 位势 z 除以 g 得到位势高度（m），与 WB2 的 Z500 定义一致。
        # Divide geopotential z by g to get geopotential height (m), as in WB2 Z500.
        obs_z500.append(regular_grid(truth_z / common.GRAVITY))
        obs_t850.append(regular_grid(truth_t))
        obs_t2m.append(regular_grid(truth_t2))
        selected_leads.append(int(leads[j]))
    leads = np.asarray(selected_leads, dtype=np.int64)
    if leads.size == 0:
        raise RuntimeError("no forecast lead has all three truth fields")

    lat = REG_LAT
    lon = REG_LON
    init = np.array([t0])
    lead_td = leads.astype("timedelta64[h]")
    forecast_ds = xr.Dataset(
        {
            # 预报场保持位势量（z），单位 m2 s-2 / keep the forecast as geopotential (m2 s-2)
            "geopotential_500": (("init_time", "prediction_timedelta", "latitude", "longitude"),
                                  np.asarray(z500)[None, :, :, :] * common.GRAVITY),
            "temperature_850": (("init_time", "prediction_timedelta", "latitude", "longitude"),
                                np.asarray(t850)[None, :, :, :]),
            "2m_temperature": (("init_time", "prediction_timedelta", "latitude", "longitude"),
                                np.asarray(t2m)[None, :, :, :]),
        },
        coords={"init_time": init, "prediction_timedelta": lead_td,
                "latitude": lat, "longitude": lon},
    )
    # 观测场以有效时刻为坐标 / observations use valid time as the coordinate
    obs_times = np.asarray([t0 + np.timedelta64(int(h), "h") for h in leads])
    obs_ds = xr.Dataset(
        {
            "geopotential_500": (("time", "latitude", "longitude"),
                                  np.asarray(obs_z500) * common.GRAVITY),
            "temperature_850": (("time", "latitude", "longitude"),
                                np.asarray(obs_t850)),
            "2m_temperature": (("time", "latitude", "longitude"), np.asarray(obs_t2m)),
        },
        coords={"time": obs_times, "latitude": lat, "longitude": lon},
    )
    forecast_ds["geopotential_500"].attrs["units"] = "m2 s-2"
    forecast_ds["temperature_850"].attrs["units"] = "K"
    forecast_ds["2m_temperature"].attrs["units"] = "K"
    obs_ds["geopotential_500"].attrs["units"] = "m2 s-2"
    obs_ds["temperature_850"].attrs["units"] = "K"
    obs_ds["2m_temperature"].attrs["units"] = "K"
    output.mkdir(parents=True, exist_ok=True)
    # WB2 tooling currently reads the conventional Zarr v2 layout.
    # 中文说明：WB2 工具目前读取常规的 Zarr v2 布局。
    forecast_ds.to_zarr(output / f"{case.key}_{method}_forecast.zarr", mode="w",
                        zarr_format=2, consolidated=True)
    obs_ds.to_zarr(output / f"{case.key}_obs.zarr", mode="w",
                   zarr_format=2, consolidated=True)
    # 同时写出 NetCDF，便于快速人工查看 / also write NetCDF for quick inspection
    forecast_ds.to_netcdf(output / f"{case.key}_{method}_forecast.nc")
    obs_ds.to_netcdf(output / f"{case.key}_obs.nc")
    print(f"wrote {output} ({leads.size} leads; bilinear 1.5° export)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-key", required=True)
    parser.add_argument("--method", choices=("direct", "nudged"), default="nudged")
    parser.add_argument("--output", type=Path, default=common.EXPERIMENT_ROOT / "data/wb2")
    args = parser.parse_args()
    case = next(case for case in common.CASES if case.key == args.case_key)
    export_case(case, args.method, args.output)


if __name__ == "__main__":
    main()
