#!/usr/bin/env python3
"""
Reproducible ERA5 -> PlaSiC T85/L25 forecast workflow.

Examples
--------
python run_workflow.py --case-key mycase --init-time 2021-03-12T00 \
    --stages download shock build run verify figures export-wb2

Model/tool binaries must be built separately with the commands shown in README.md.
"""

from __future__ import annotations
import argparse
import datetime as dt
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

SCRIPTS = ROOT / "scripts"

sys.path.insert(0, str(SCRIPTS))


# ----------------------------------------------------------------------
# 导入工作流各阶段对应的本地模块
# Import local modules corresponding to individual workflow stages
# ----------------------------------------------------------------------

# Build PlaSiC initial conditions.
import build_ic  # noqa: E402

# Contains shared configuration, the ForecastCase definition, archived cases,
# and nudging-related parameters.
import common  # noqa: E402

# Download ERA5 data from the Copernicus Climate Data Store.
import download_era5  # noqa: E402

# Download shock-related data.
import download_shock  # noqa: E402

# Extract the numeric report tables from the verification archives.
import make_report_tables  # noqa: E402

# Run PlaSiC forecast integrations.
import run_forecasts  # noqa: E402

# Verify model forecast results.
import verify  # noqa: E402

# Export forecast output to a WeatherBench2-compatible format.
import export_wb2  # noqa: E402


def selected_case(args: argparse.Namespace) -> common.ForecastCase:
    """
    根据命令行参数确定本次工作流需要处理的预报个例。
    Determine which forecast case should be processed from command-line arguments.

    参数
    Parameters
    ----------
    args : argparse.Namespace
        argparse 解析得到的命令行参数。
        Command-line arguments parsed by argparse.

    返回
    Returns
    -------
    common.ForecastCase
        A ForecastCase object containing the case key, initialization time,
        label, description, and related metadata.

    逻辑
    Logic
    -----
    1. 如果没有提供 --init-time，则从 common.CASES 中选择归档个例。
       If --init-time is absent, select an archived case from common.CASES.

    2. 如果提供了 --init-time，则动态创建一个新的 ForecastCase。
       If --init-time is provided, dynamically construct a new ForecastCase.
    """

    # ------------------------------------------------------------------
    # 情况 1：用户没有提供新的初始化时间
    # Case 1: no new initialization time was supplied
    # ------------------------------------------------------------------
    if not args.init_time:

        # 从归档个例 common.CASES 中筛选候选项。
        # Filter candidate cases from the archived common.CASES list.
        matches = [
            c for c in common.CASES if not args.case_key or c.key == args.case_key
        ]

        # 必须最终唯一确定一个个例。
        # Exactly one case must be identified.
        if len(matches) != 1:
            raise SystemExit(
                "Specify --case-key for an archived case "
                "or provide --init-time for a new case"
            )

        return matches[0]

    # ------------------------------------------------------------------
    # 情况 2：用户通过 --init-time 指定了一个新的个例
    # Case 2: a new case was specified via --init-time
    # ------------------------------------------------------------------

    # 将 ISO 格式字符串解析为 datetime 对象。
    # Parse the ISO-formatted initialization time into a datetime object.
    #
    # 示例 / Example:
    #   "2021-03-12T00"
    #       ->
    #   datetime.datetime(2021, 3, 12, 0, 0)
    when = dt.datetime.fromisoformat(args.init_time)

    # 如果用户指定了 --case-key，则直接使用它。
    # Use the user-specified --case-key when available.
    #
    # 否则根据初始化时刻自动生成 case key，例如：
    # Otherwise generate one from the initialization time, e.g.:
    #
    #   case_2021031200
    key = args.case_key or when.strftime("case_%Y%m%d%H")

    # 为新的初始化时刻创建 ForecastCase。
    # Construct a ForecastCase for the new initialization time.
    return common.ForecastCase(
        # 个例唯一标识符。
        # Unique identifier of the forecast case.
        key=key,
        # 统一保存为 YYYY-MM-DDTHH 格式。
        # Store the initialization time consistently as YYYY-MM-DDTHH.
        init_time=when.strftime("%Y-%m-%dT%H"),
        # 个例显示标签。
        # Human-readable case label.
        label_zh=f"{key}: {when:%Y-%m-%d %H:%M} UTC",
        # 对用户自定义个例进行简单说明。
        # Brief description for a user-selected initialization case.
        description_zh="User-selected ERA5 initialisation case.",
    )


def run_case(
    case: common.ForecastCase,
    stages: set[str],
    force: bool,
    nudge_window_hours: int | None = None,
    nudge_tau_hours: float | None = None,
    surface_tau_hours: float | None = None,
    require_exact: bool = False,
) -> None:
    """
    按指定阶段执行单个预报个例的完整工作流。
    Execute selected workflow stages for a single forecast case.

    参数
    Parameters
    ----------
    case : common.ForecastCase
        当前需要运行的预报个例。
        Forecast case to process.

    stages : set[str]
        要执行的工作流阶段集合，例如：
        Set of workflow stages to execute

    force : bool
        是否强制重新执行已经存在的模型积分结果。
        Whether to force rerunning model segments that already exist.

    nudge_window_hours : int | None
        ERA5 nudging/relaxation 窗口长度，单位为小时。
        Length of the ERA5 nudging/relaxation window in hours.

        如果为 None，则沿用 common 模块中的默认值。
        If None, retain the default configured in common.

    nudge_tau_hours : float | None
        大气变量 Newtonian relaxation 的时间尺度，单位为小时。
        Atmospheric Newtonian-relaxation time scale in hours.

        如果为 None，则使用默认值。
        If None, retain the configured default.

    surface_tau_hours : float | None
        地表变量 relaxation 时间尺度，单位为小时。
        Surface-variable relaxation time scale in hours.

        如果为 None，则使用默认值。
        If None, retain the configured default.

    require_exact : bool
        检验阶段是否严格要求 ERA5 时间轴完整。
        Whether verification should require complete ERA5 time axes.
    """

    # ------------------------------------------------------------------
    # 限制当前运行仅处理传入的 case
    # Restrict the workflow to the supplied case
    # ------------------------------------------------------------------

    common.CASES = [case]

    # ------------------------------------------------------------------
    # 可选：覆盖 nudging 窗口长度
    # Optional: override the nudging-window duration
    # ------------------------------------------------------------------
    if nudge_window_hours is not None:

        # nudging 窗口必须至少为 1 小时。
        # The nudging window must be at least one hour.
        if nudge_window_hours < 1:
            raise ValueError("--nudge-window-hours must be positive")

        common.NUDGE_WINDOW_HOURS = nudge_window_hours

    # ------------------------------------------------------------------
    # 可选：覆盖大气 nudging 时间尺度
    # Optional: override the atmospheric nudging time scale
    # ------------------------------------------------------------------
    if nudge_tau_hours is not None:

        # relaxation 时间尺度必须严格大于 0。
        # The relaxation time scale must be strictly positive.
        if nudge_tau_hours <= 0:
            raise ValueError("--tau-hours must be positive")

        common.NUDGE_TAU_HOURS = nudge_tau_hours

    # ------------------------------------------------------------------
    # 可选：覆盖地表 nudging 时间尺度
    # Optional: override the surface nudging time scale
    # ------------------------------------------------------------------
    if surface_tau_hours is not None:
        if surface_tau_hours <= 0:
            raise ValueError("--surface-tau-hours must be positive")

        common.NUDGE_SURFACE_TAU_HOURS = surface_tau_hours

    # ------------------------------------------------------------------
    # Stage: download
    # 阶段：下载 ERA5
    # ------------------------------------------------------------------
    if "download" in stages:

        import cdsapi

        client = cdsapi.Client(
            retry_max=5,
            sleep_max=30,
        )

        download_era5.run_case(client, case)

    # ------------------------------------------------------------------
    # Stage: shock
    # 阶段：下载 shock
    # ------------------------------------------------------------------
    if "shock" in stages:

        # 下载或准备当前 case 对应的 shock 数据。
        # Download or prepare shock-related data for the current case.
        download_shock.run_case(case)

    # ------------------------------------------------------------------
    # Stage: build
    # 阶段：构建 PlaSiC 初始条件
    # ------------------------------------------------------------------
    if "build" in stages:

        build_ic.build_case(
            case,
            nudge=True,
        )

    # ------------------------------------------------------------------
    # Stage: run
    # 阶段：执行 PlaSiC 模型预报
    # ------------------------------------------------------------------
    if "run" in stages:

        run_forecasts.run_case(
            case,
            ["direct", "nudged"],
            force=force,
        )

    # ------------------------------------------------------------------
    # Stage: verify
    # 阶段：预报检验
    # ------------------------------------------------------------------
    if "verify" in stages:

        verify.VERIF_ROOT.mkdir(
            parents=True,
            exist_ok=True,
        )

        import numpy as np

        for method in ("direct", "nudged"):

            result = verify.verify_case(
                case,
                method,
                require_exact=require_exact,
            )
            np.savez_compressed(
                verify.VERIF_ROOT / f"{case.key}_{method}.npz",
                **result,
            )

    # ------------------------------------------------------------------
    # Stage: figures
    # 阶段：生成图像
    # ------------------------------------------------------------------
    if "figures" in stages:

        import make_figures
        import make_nature_figures

        make_nature_figures.CASES = [case.key]

        make_nature_figures.CASE_LABELS[case.key] = case.key

        make_nature_figures.main()

        make_figures.main_for_cases([case])

    # ------------------------------------------------------------------
    # Stage: export-wb2
    # 阶段：导出 WeatherBench2 数据
    # ------------------------------------------------------------------
    if "export-wb2" in stages:

        export_wb2.export_case(
            case,
            "direct",
            ROOT / "data/wb2",
        )
        export_wb2.export_case(
            case,
            "nudged",
            ROOT / "data/wb2",
        )


def main() -> None:
    """
    程序命令行入口。
    Command-line entry point.

    该函数负责：
    This function is responsible for:

    1. 定义命令行参数；
       Defining command-line arguments;

    2. 解析用户输入；
       Parsing user input;

    3. 确定当前预报个例；
       Determining the forecast case;

    4. 调用 run_case() 执行所选工作流阶段。
       Calling run_case() to execute the selected workflow stages.
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--case-key",
        help="archived case key, or name for a new case",
    )

    # ------------------------------------------------------------------
    # --init-time
    # ------------------------------------------------------------------

    parser.add_argument(
        "--init-time",
        help="new case ISO time, e.g. 2021-03-12T00",
    )

    # ------------------------------------------------------------------
    # --stages
    # ------------------------------------------------------------------

    parser.add_argument(
        "--stages",
        nargs="+",
        default=[
            "build",
            "run",
            "verify",
            "figures",
        ],
        choices=[
            "download",
            "shock",
            "build",
            "run",
            "verify",
            "figures",
            "export-wb2",
        ],
    )

    # ------------------------------------------------------------------
    # --force
    # ------------------------------------------------------------------

    parser.add_argument(
        "--force",
        action="store_true",
        help="rerun existing model segments",
    )

    # ------------------------------------------------------------------
    # --nudge-window-hours
    # ------------------------------------------------------------------

    parser.add_argument(
        "--nudge-window-hours",
        type=int,
        default=None,
        help="ERA5 relaxation window (default: 12 h)",
    )

    # ------------------------------------------------------------------
    # --tau-hours
    # ------------------------------------------------------------------

    parser.add_argument(
        "--tau-hours",
        type=float,
        default=None,
        help="atmospheric Newtonian relaxation time scale (default: 1 h)",
    )

    # ------------------------------------------------------------------
    # --surface-tau-hours
    # ------------------------------------------------------------------

    parser.add_argument(
        "--surface-tau-hours",
        type=float,
        default=None,
        help="surface relaxation time scale (default: 24 h)",
    )

    # ------------------------------------------------------------------
    # --require-exact
    # ------------------------------------------------------------------

    parser.add_argument(
        "--require-exact",
        action="store_true",
        help=("require complete CDS hourly and six-hourly ERA5 axes " "during verify"),
    )

    # ------------------------------------------------------------------
    # 解析命令行
    # Parse command-line arguments
    # ------------------------------------------------------------------

    args = parser.parse_args()

    run_case(
        selected_case(args),
        set(args.stages),
        args.force,
        args.nudge_window_hours,
        args.tau_hours,
        args.surface_tau_hours,
        args.require_exact,
    )


if __name__ == "__main__":
    main()
