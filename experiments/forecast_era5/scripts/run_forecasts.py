#!/usr/bin/env python3
"""Run the PlaSiC 10-day forecasts (direct and nudged initial conditions).

中文说明：运行 PlaSiC 10 天预报（直接插值初值与松弛同化初值两种方案）。

Every forecast is integrated in two segments so that the frame stream carries
both an hourly view of the first day (for the initial-shock analysis) and the
6-hourly view of the remaining nine days:

中文说明：每次预报分两段积分，使帧流既包含首日逐小时视图（用于初始冲击
分析），又包含余下九天每 6 小时的视图：

* segment 1: 0 -> 24 h, one frame per hour (``--frame-interval 4``);
  中文说明：第 1 段：0 → 24 h，每小时输出一帧（``--frame-interval 4``）；
* segment 2: 24 h -> 240 h, one frame per 6 h (``--frame-interval 24``).
  中文说明：第 2 段：24 h → 240 h，每 6 小时输出一帧（``--frame-interval 24``）。

Two initial conditions are compared for every case:

中文说明：每个个例比较两种初始条件：

* ``direct``  the ERA5 analysis at t0, mapped to T85/L25 (old time level equal
  to the current one);
  中文说明：``direct``——t0 时刻的 ERA5 分析场映射到 T85/L25（蛙跳格式的
  旧时间层等于当前时间层）；
* ``nudged``  a 12-h Newtonian-relaxation initialization (tau = 1 h) toward the
  hourly ERA5 analyses, implemented by 1-h restart recycling, ending at t0.
  中文说明：``nudged``——向逐小时 ERA5 分析场做 12 小时牛顿松弛初始化
  （τ = 1 h），通过每小时循环一次重启文件实现，终点为 t0。
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_ic  # noqa: E402
import common  # noqa: E402
import frame_io  # noqa: E402
import restart_io  # noqa: E402

# 模式可执行文件路径可由环境变量覆盖，默认使用 T85/L25 MPI P4 构建。
# The model executable can be overridden by an environment variable; the
# default is the T85/L25 MPI P4 build.
PROGRAM = Path(os.environ.get(
    "PLASIC_EXECUTABLE",
    str(common.PROJECT_ROOT / "build/T85_L25_MPI_P4_O16/plasic.x"),
))
NPRO = int(os.environ.get("PLASIC_NPRO", "4"))  # MPI 进程数 / number of MPI ranks
SEGMENT1_STEPS = common.NTSPD                      # 第一段 24 小时 / segment 1: 24 h
SEGMENT1_FRAME_INTERVAL = common.STEPS_PER_HOUR    # 第一段帧间隔 1 小时 / segment 1 frame interval: 1 h
SEGMENT2_STEPS = common.FORECAST_STEPS - SEGMENT1_STEPS  # 第二段 216 小时 / segment 2: remaining hours
SEGMENT2_FRAME_INTERVAL = common.STEPS_PER_6H      # 第二段帧间隔 6 小时 / segment 2 frame interval: 6 h


def run_model(restart: Path, steps: int, output: Path, frames: Path | None,
              frame_interval: int, log: Path, co2_ppm: float = 410.0) -> None:
    """Run one model segment from a restart file.

    中文说明：从重启文件出发运行一段模式积分。

    ``frames`` may be ``None`` for the short nudging cycles, which do not need
    a frame stream; stdout/stderr are redirected to ``*.stdout`` next to the
    progress log.
    中文说明：松弛循环等短积分不需要帧流，此时 ``frames`` 可为 ``None``；
    标准输出/标准错误重定向到进度日志同名的 ``*.stdout`` 文件。
    """
    log.parent.mkdir(parents=True, exist_ok=True)
    # 单进程直接运行，多进程通过 mpirun 启动 / direct run for one rank, mpirun otherwise
    launcher = [str(PROGRAM)] if NPRO <= 1 else ["mpirun", "-np", str(NPRO), str(PROGRAM)]
    command = launcher + [
        "--restart", str(restart),
        "--co2-ppm", repr(co2_ppm),
        "--steps", str(steps),
        "--output", str(output),
        "--progress", str(log),
        "--progress-interval", str(max(1, steps // 8)),
    ]
    if frames is not None:
        command += ["--frames", str(frames), "--frame-interval", str(frame_interval)]
    started = time.time()
    with log.with_suffix(".stdout").open("w") as stream:
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)
    print(f"    {steps} steps in {time.time() - started:.0f} s -> {output.name}")


def blend_restart(model_restart: Path, target: dict[str, np.ndarray],
                  gamma: float, output: Path) -> None:
    """Insert the analysis increment into the current time level only.

    中文说明：只把分析增量写入当前时间层。

    The model's own previous (filtered) leapfrog level is kept, so the
    integration stays continuous and the Robert-Asselin filter damps the
    increment smoothly instead of being restarted from a cold state every hour.

    中文说明：保留模式自身的上一（已滤波）蛙跳时间层，使积分保持连续，
    Robert–Asselin 滤波器能平滑地抑制增量，而不是每小时都从“冷启动”重新开始。
    """
    records = restart_io.read_raw_records(model_restart)
    replacements: dict[str, bytes] = {}
    for name in ("st", "sd", "sz", "sq", "sp"):
        model = restart_io.decode_float(dict(records)[name])
        # 新值 = 目标 + γ·(模式 − 目标)，γ = exp(−Δt/τ) 由调用方计算。
        # new = target + gamma * (model - target), with gamma = exp(-dt/tau).
        blended = np.asarray(target[name], dtype="<f4") + gamma * (
            model - np.asarray(target[name], dtype="<f4")
        )
        replacements[name] = restart_io.encode_float(blended)
    restart_io.write_raw_records(
        output, restart_io.replace_records(records, replacements)
    )


def blend_surface_restart(model_restart: Path, target: dict[str, np.ndarray],
                          gamma: float, output: Path) -> None:
    """Relax slowly varying surface records toward the matching ERA5 hour.

    中文说明：把缓变的地表状态记录向对应时刻的 ERA5 场松弛。

    Missing/land sentinel values are copied from the model, so a finite
    analysis never contaminates the model's land/ocean bookkeeping.

    中文说明：缺测值/陆地哨兵值直接从模式状态复制，确保有限的分析值不会污染
    模式对陆地/海洋状态的记账。
    """
    records = restart_io.read_raw_records(model_restart)
    current = dict(records)
    replacements: dict[str, bytes] = {}
    for name, target_values in target.items():
        if name not in current:
            continue
        model = restart_io.decode_float(current[name]).astype(np.float64)
        target_array = np.asarray(target_values, dtype=np.float64).reshape(model.shape)
        # 哨兵值（≥1e10 或非有限值）不参与松弛 / sentinels (>=1e10 or non-finite) are not relaxed
        valid = np.isfinite(target_array) & (np.abs(target_array) < 1.0e10)
        blended = np.where(valid, target_array + gamma * (model - target_array), model)
        replacements[name] = restart_io.encode_float(blended)
    restart_io.write_raw_records(output, restart_io.replace_records(records, replacements))


def frames_match(path: Path, expected_count: int, expected_interval: int) -> bool:
    """Return true only for a frame stream with the requested cadence.

    中文说明：仅当帧流具有所要求的输出节奏时才返回 True。

    This prevents a pre-6-hour archive from being silently reused after the
    output schedule changes.
    中文说明：这样可以防止旧的（非 6 小时节奏的）归档在输出计划变更后被
    悄悄复用。
    """
    if not path.exists():
        return False
    try:
        grouped = frame_io.load_stream(path)
        frames = next(iter(grouped.values()))
        steps = np.asarray([frame.step for frame in frames], dtype=np.int64)
    except (OSError, ValueError, StopIteration):
        return False
    return (steps.size == expected_count
            and (steps.size < 2 or np.all(np.diff(steps) == expected_interval)))


def forecast(ic_path: Path, run_dir: Path, frame_dir: Path, tag: str,
             force: bool) -> None:
    """Run the two-segment 10-day forecast for one initial condition.

    中文说明：对一种初始条件运行两段式 10 天预报。
    """
    frames1 = frame_dir / f"{tag}_seg1.frames"
    frames2 = frame_dir / f"{tag}_seg2.frames"
    # 预期帧数 = 积分步数 ÷ 帧间隔 / expected frames = steps / frame interval
    expected1 = SEGMENT1_STEPS // SEGMENT1_FRAME_INTERVAL
    expected2 = SEGMENT2_STEPS // SEGMENT2_FRAME_INTERVAL
    if (frames_match(frames1, expected1, SEGMENT1_FRAME_INTERVAL)
            and frames_match(frames2, expected2, SEGMENT2_FRAME_INTERVAL)
            and not force):
        print(f"    [skip] {tag} frames exist")
        return
    print(f"  [{tag}] segment 1: 0-24 h, hourly frames")
    run_model(ic_path, SEGMENT1_STEPS, run_dir / f"{tag}_seg1.restart",
              frames1, SEGMENT1_FRAME_INTERVAL, run_dir / f"{tag}_seg1.jsonl")
    print(f"  [{tag}] segment 2: 24-240 h, 6-hourly frames")
    # 第二段从第一段结束时的重启文件继续。/ segment 2 continues from segment 1's restart.
    run_model(run_dir / f"{tag}_seg1.restart", SEGMENT2_STEPS,
              run_dir / f"{tag}_final.restart", frames2, SEGMENT2_FRAME_INTERVAL,
              run_dir / f"{tag}_seg2.jsonl")


def nudged_ic(case: common.ForecastCase, run_dir: Path, ic_dir: Path, force: bool = False) -> Path:
    """Build the nudged initial condition by 1-hour restart recycling.

    中文说明：通过每小时循环重启文件构建松弛同化初始条件。

    Each hour: run 1 h from the current state, relax the atmospheric spectra
    toward the hourly ERA5 target with weight gamma, then relax the surface
    records with the slower surface weight.  The final state is written as
    ``nudged.restart``.
    中文说明：每小时先积分 1 小时，再按权重 gamma 把大气谱系数向对应时刻的
    ERA5 目标松弛，然后用更慢的地表权重松弛地表记录。最终状态写为
    ``nudged.restart``。
    """
    nudged_ic_path = ic_dir / "nudged.restart"
    if nudged_ic_path.exists() and not force:
        print(f"    [skip] {nudged_ic_path.name} exists")
        return nudged_ic_path
    print(f"  nudging initialization ({common.NUDGE_WINDOW_HOURS} h, "
          f"tau = {common.NUDGE_TAU_HOURS} h)")
    entries = build_ic.read_nudge_file(ic_dir / "nudge_targets.bin")
    surface_path = ic_dir / "surface_targets.npz"
    surface_entries = build_ic.read_surface_targets(surface_path) if surface_path.exists() else []
    # 大气谱的松弛权重：γ = exp(−Δt/τ)，Δt = 1 h / atmospheric relaxation weight
    gamma = float(np.exp(-1.0 / common.NUDGE_TAU_HOURS))
    state = ic_dir / "nudge_start.restart"
    for hour in range(1, len(entries)):
        step, target = entries[hour]
        cycle_out = run_dir / f"nudge_cycle_{hour:02d}.restart"
        # 从当前状态向前积分 1 小时 / integrate one hour from the current state
        run_model(state, common.STEPS_PER_HOUR, cycle_out, None, 1,
                  run_dir / f"nudge_cycle_{hour:02d}.jsonl")
        blended = run_dir / f"nudge_blend_{hour:02d}.restart"
        blend_restart(cycle_out, target, gamma, blended)
        if surface_entries:
            surface_target = surface_entries[hour][1]
            surface_blended = run_dir / f"nudge_surface_blend_{hour:02d}.restart"
            # 地表状态使用更长的松弛时间尺度（24 h）。/ slower surface relaxation (24 h).
            surface_gamma = float(np.exp(-1.0 / common.NUDGE_SURFACE_TAU_HOURS))
            blend_surface_restart(blended, surface_target, surface_gamma, surface_blended)
            state = surface_blended
        else:
            state = blended
        # 校验循环结束时的模式步数确实等于目标时刻，防止时间轴错位。
        # Verify the cycle ends exactly at the target step to catch clock drift.
        nstep = restart_io.decode_int(
            dict(restart_io.read_raw_records(state))["nstep"]
        )
        assert nstep == step, (nstep, step)
    restart_io.write_raw_records(nudged_ic_path, restart_io.read_raw_records(state))
    print(f"    nudged IC written -> {nudged_ic_path.name}")
    return nudged_ic_path


def run_case(case: common.ForecastCase, methods: "list[str]", force: bool) -> None:
    """Run the requested initialization methods for one case.

    中文说明：对单个个例运行所要求的初值方案。
    """
    run_dir = common.RUN_ROOT / case.key
    frame_dir = common.FRAME_ROOT / case.key
    ic_dir = common.IC_ROOT / case.key
    run_dir.mkdir(parents=True, exist_ok=True)
    frame_dir.mkdir(parents=True, exist_ok=True)

    if "direct" in methods:
        print(f"[{case.key}] direct forecast")
        forecast(ic_dir / "direct.restart", run_dir, frame_dir, "direct", force)
    if "nudged" in methods:
        print(f"[{case.key}] nudged forecast")
        ic = nudged_ic(case, run_dir, ic_dir, force=force)
        forecast(ic, run_dir, frame_dir, "nudged", force)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--methods", nargs="*", default=["direct", "nudged"])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        run_case(case, args.methods, args.force)


if __name__ == "__main__":
    main()
