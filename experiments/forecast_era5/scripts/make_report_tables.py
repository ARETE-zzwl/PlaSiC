#!/usr/bin/env python3
"""Extract the numeric tables in the report from the verification archives.

中文说明：从检验归档文件中提取报告所需的数值表格（以 Markdown 打印）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

VERIF_ROOT = common.DATA_ROOT / "verification"  # 检验归档目录 / verification archive directory


def load(case_key: str, method: str) -> dict:
    """Load one case/method verification archive.

    中文说明：载入某个“个例 + 初值方法”的检验归档文件。
    """
    return dict(np.load(VERIF_ROOT / f"{case_key}_{method}.npz", allow_pickle=True))


def at_lead(data: dict, key: str, lead_hour: float, level_index=None, scale=1.0):
    """Pick the archived value nearest to ``lead_hour``.

    中文说明：取归档数据中最接近 ``lead_hour`` 的预报时效对应的值。

    ``level_index`` selects a level if the metric is level-dependent and
    ``scale`` converts units.
    中文说明：``level_index`` 用于选取层次（当指标与层次有关时），``scale``
    用于单位换算。
    """
    index = int(np.argmin(np.abs(data["all_lead_hours"] - lead_hour)))
    values = data[key]
    values = values[level_index] if level_index is not None else values
    return values[index] * scale


def table_rmse() -> None:
    """Print the per-case/per-method RMSE summary table (Markdown).

    中文说明：打印按个例、初值方法、预报时效分列的 RMSE 汇总表（Markdown）。
    """
    # 中文说明：每项为 (表头, 归档键, 层次索引, 单位换算, 单位)。
    # English: each spec is (label, archive key, level index, scale, unit).
    specs = [
        ("T850", "rmse_t", 0, 1.0, "K"),
        ("T500", "rmse_t", 1, 1.0, "K"),
        ("T200", "rmse_t", 2, 1.0, "K"),
        ("u250", "rmse_u", 1, 1.0, "m/s"),
        ("q700", "rmse_q", 1, 1000.0, "g/kg"),
        ("Z500", "rmse_z", 0, 1.0, "m"),
        ("T2m", "rmse_t2m", None, 1.0, "K"),
        ("MSLP", "rmse_mslp", None, 0.01, "hPa"),
    ]
    leads = [24, 72, 120, 240]
    print("### RMSE 汇总（每个个例、每种初值）\n")
    header = "| 变量 | " + " | ".join(
        f"{case.key}-{method}-{lead}h" for case in common.CASES
        for method in ("direct", "nudged") for lead in leads
    ) + " |"
    print(header)
    print("|" + "---|" * (len(header.split("|")) - 2))
    for label, key, level, scale, unit in specs:
        row = []
        for case in common.CASES:
            for method in ("direct", "nudged"):
                data = load(case.key, method)
                for lead in leads:
                    row.append(f"{at_lead(data, key, lead, level, scale):.2f}")
        print(f"| {label} ({unit}) | " + " | ".join(row) + " |")


def table_case_mean() -> None:
    """Print the three-case-mean RMSE table (Markdown).

    中文说明：打印三个个例平均后的 RMSE 表（Markdown）。
    """
    # 中文说明：每项为 (表头, 归档键, 层次索引, 单位换算)。
    # English: each spec is (label, archive key, level index, scale).
    specs = [
        ("T850 (K)", "rmse_t", 0, 1.0),
        ("T500 (K)", "rmse_t", 1, 1.0),
        ("u250 (m/s)", "rmse_u", 1, 1.0),
        ("q700 (g/kg)", "rmse_q", 1, 1000.0),
        ("Z500 (m)", "rmse_z", 0, 1.0),
        ("MSLP (hPa)", "rmse_mslp", None, 0.01),
    ]
    leads = [24, 48, 72, 120, 168, 240]
    print("\n### 三个个例平均 RMSE\n")
    header = "| 变量 | " + " | ".join(
        f"{method}-{lead}h" for method in ("direct", "nudged") for lead in leads
    ) + " |"
    print(header)
    print("|" + "---|" * (len(header.split("|")) - 2))
    for label, key, level, scale in specs:
        row = []
        for method in ("direct", "nudged"):
            for lead in leads:
                values = []
                for case in common.CASES:
                    data = load(case.key, method)
                    values.append(at_lead(data, key, lead, level, scale))
                # 对不同个例取平均（自动忽略 NaN）/ average over cases, ignoring NaN
                row.append(f"{np.nanmean(values):.2f}")
        print(f"| {label} | " + " | ".join(row) + " |")


def table_shock() -> None:
    """Print the initial-shock / initial-condition-quality table.

    中文说明：打印“初始冲击 / 初值质量”表（以 2020-01-15 个例为例）。
    """
    print("\n### 初始冲击与初值质量（个例 A：2020-01-15）\n")
    print("| 指标 | 直接插值 | 松弛同化 |")
    print("|---|---:|---:|")
    rows = []
    for method in ("direct", "nudged"):
        data = load("winter2020", method)
        # 第 0 个时效（t0 时刻）的 RMSE：检验初值本身的偏差。
        # RMSE at lead index 0 (t0): measures the analysis increment itself.
        rows.append(
            (
                f"{data['rmse_t'][1][0]:.3f}",
                f"{data['rmse_t'][0][0]:.3f}",
                f"{data['rmse_u'][1][0]:.3f}",
                f"{data['rmse_q'][1][0] * 1000:.4f}",
            )
        )
    print(f"| t₀ 时刻 T500 与 ERA5 的 RMSE (K) | {rows[0][0]} | {rows[1][0]} |")
    print(f"| t₀ 时刻 T850 与 ERA5 的 RMSE (K) | {rows[0][1]} | {rows[1][1]} |")
    print(f"| t₀ 时刻 u250 与 ERA5 的 RMSE (m/s) | {rows[0][2]} | {rows[1][2]} |")
    print(f"| t₀ 时刻 q700 与 ERA5 的 RMSE (g/kg) | {rows[0][3]} | {rows[1][3]} |")
    for method in ("direct", "nudged"):
        data = load("winter2020", method)
        # 初值质量在起报后 1/6/12/24 小时的衰减。
        # Decay of the initial-condition quality at +1/6/12/24 h.
        print(
            f"| {method}: T500 RMSE 1 h / 6 h / 12 h / 24 h (K) | "
            f"{at_lead(data, 'rmse_t', 1, 1):.2f} / {at_lead(data, 'rmse_t', 6, 1):.2f} / "
            f"{at_lead(data, 'rmse_t', 12, 1):.2f} / {at_lead(data, 'rmse_t', 24, 1):.2f} |"
        )


def table_persistence() -> None:
    """Print the persistence-baseline comparison table.

    中文说明：打印与持续性（persistence）预报基线对比的表格。
    """
    print("\n### 与持续性预报对比（三个个例平均）\n")
    print("| 变量 | 预报时效 | 直接插值 | 持续性 |")
    print("|---|---:|---:|---:|")
    for label, key, level, scale in [
        ("T500 (K)", "rmse_t", 1, 1.0),
        ("Z500 (m)", "rmse_z", 0, 1.0),
        ("u250 (m/s)", "rmse_u", 1, 1.0),
    ]:
        for lead in (48, 120, 240):
            model, persist = [], []
            for case in common.CASES:
                data = load(case.key, "direct")
                model.append(at_lead(data, key, lead, level, scale))
                # 持续性基线在归档中的键名为 persistence_t500 / persistence_z500 / persistence_u250。
                # Persistence baselines are archived as
                # persistence_t500 / persistence_z500 / persistence_u250.
                persist_key = f"persistence_{'t500' if key == 'rmse_t' else ('z500' if key == 'rmse_z' else 'u250')}"
                if persist_key in data:
                    persist.append(at_lead(data, persist_key, lead, None, scale))
            persist_value = np.nanmean(persist) if persist else float("nan")
            print(f"| {label} | {lead} h | {np.nanmean(model):.2f} | {persist_value:.2f} |")


def main() -> None:
    table_rmse()
    table_case_mean()
    table_shock()
    table_persistence()


if __name__ == "__main__":
    main()
