#!/usr/bin/env python3
"""Extra hourly ERA5 fields for the initial-shock analysis of one case.

中文说明：为单个个例的“初始冲击（initial shock）”分析补充逐小时 ERA5 场。

English: downloads t/u/v/q/z on the IC pressure levels plus the single-level
fields for t0 .. t0+24 h (hourly), so every archived case has a complete
first-day verification stream.

中文说明：下载 t0 至 t0+24 h 逐小时的、IC 气压层上的 t/u/v/q/z 以及单层
变量场，使每个已归档个例都有完整的首日检验序列。
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import download_era5 as base  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None,
                        help="case keys (default: all)")
    parser.add_argument("--force", action="store_true",
                        help="replace existing shock files, for example to add geopotential")
    args = parser.parse_args()
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        run_case(case, force=args.force)


def run_case(case, force: bool = False) -> None:
    # cdsapi 只在真正下载时才导入，使 build/run/verify 不依赖它。
    # Import cdsapi only when a download actually runs, so build/run/verify
    # do not depend on it.
    import cdsapi

    t0 = common.init_datetime(case)
    directory = common.ERA5_ROOT / case.key
    directory.mkdir(parents=True, exist_ok=True)
    # 初始冲击分析只关心前 25 个整点（t0 到 t0+24 h）。
    # The initial-shock analysis only needs the first 25 exact hours (t0..t0+24 h).
    times = [t0 + dt.timedelta(hours=h) for h in range(0, 25)]
    client = cdsapi.Client(retry_max=5, sleep_max=30)
    # Older archives contain only t/u/v/q.  Re-fetch them once so hourly
    # geopotential-height verification is available too.
    # 中文说明：旧归档只含 t/u/v/q，这里自动检测并重新下载一次，以便支持
    # 逐小时位势高度检验。
    if not force:
        import xarray as xr
        for existing in sorted(directory.glob("shock_pl_*.nc")):
            try:
                with xr.open_dataset(existing) as dataset:
                    if "geopotential" not in dataset.data_vars and "z" not in dataset.data_vars:
                        force = True
                        break
            except (OSError, ValueError):
                force = True  # 无法读取的文件也视为需要重下 / unreadable files are re-downloaded too

    # 气压层场：按“天”分组请求，避免 CDS 的 day×time 笛卡尔积歧义。
    # Pressure-level fields: request one day at a time to avoid the CDS
    # day-by-time Cartesian product silently adding unrequested hours.
    for day in sorted({when.date() for when in times}):
        hours = [when for when in times if when.date() == day]
        request = {
            "product_type": ["reanalysis"],
            "variable": base.PL_VARIABLES_VERIF,
            "pressure_level": [str(level) for level in common.ERA5_IC_LEVELS],
            "data_format": "netcdf",
            "download_format": "unarchived",
            "area": base.AREA,
            "grid": [common.ERA5_GRID, common.ERA5_GRID],
            **base.request_times(hours),
        }
        base.download(client, base.DATASET_PL, request,
                      directory / f"shock_pl_{day:%Y%m%d}.nc", force=force,
                      expected_times=hours,
                      required_variables=base.PL_VARIABLES_VERIF,
                      pressure_levels=common.ERA5_IC_LEVELS)

    # Keep the request daily.  CDS interprets ``day`` and ``time`` as a
    # Cartesian product, so one request spanning two dates would silently
    # return 48 timestamps instead of the requested 25.  The verifier reads
    # the resulting ``shock_surface_YYYYMMDD.nc`` files as one time series.
    # 中文说明：请求必须按天拆分。CDS 把 ``day`` 和 ``time`` 解释为笛卡尔
    # 积，跨两天的单个请求会悄悄返回 48 个时刻而不是请求的 25 个。检验脚本
    # 会把生成的 ``shock_surface_YYYYMMDD.nc`` 文件当作一条时间序列读取。
    base_request = {
        "product_type": ["reanalysis"],
        # Shock verification needs 2-m temperature, surface pressure and
        # mean-sea-level pressure.  The IC-only list omits MSLP, so use the
        # verification list here rather than silently reconstructing it.
        # 中文说明：冲击检验需要 2 m 温度、地面气压和海平面气压。仅含初值的
        # 变量列表没有 MSLP，因此这里直接使用检验变量列表，而不是事后重构。
        "variable": base.SL_VARIABLES_VERIF,
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": base.AREA,
        "grid": [common.ERA5_GRID, common.ERA5_GRID],
    }
    base._download_exact_days(
        client,
        base.DATASET_SL,
        base_request,
        times,
        directory / "shock_surface.nc",
        force=force,
        required_variables=base.SL_VARIABLES_VERIF,
    )


if __name__ == "__main__":
    main()
