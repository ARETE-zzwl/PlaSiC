#!/usr/bin/env python3
"""Train a tiny next-step emulator on PlaSiC output.

中文说明：在 PlaSiC 输出上训练小型 next-step emulator 的示例脚本。
包含两个基线（持续性预报、气候态）和一个小型 CNN；用 mask 掉掩膜点后的
MSE 评估。先在合成数据上跑通，再换成真实的 PlaSiC 输出。

English: baselines (persistence, climatology) + a small CNN emulator,
train/val split by time, masked MSE. Run on synthetic data first, then on
real PlaSiC output.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from plasic_dataset import PlasicDataset


class TinyCNN(nn.Module):
    """Small fully-convolutional next-step model. 中文说明：小型全卷积单步模式。"""

    def __init__(self, channels: int, hidden: int = 32) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(channels, hidden, 3, padding=1), nn.ReLU(),
            nn.Conv2d(hidden, hidden, 3, padding=1), nn.ReLU(),
            nn.Conv2d(hidden, channels, 3, padding=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)  # residual: predict the increment


def masked_mse(pred: torch.Tensor, target: torch.Tensor,
               mask: torch.Tensor) -> torch.Tensor:
    """MSE over valid grid points only. 中文说明：只在有效格点上算 MSE。"""
    se = (pred - target) ** 2
    return (se * mask).sum() / mask.sum().clamp(min=1)


@torch.no_grad()
def evaluate(baseline: str, loader: DataLoader) -> float:
    """Baseline MSE: 'persistence' (y_hat = x) or 'climatology' (y_hat = mean).

    中文说明：基线评估，persistence 用输入当预测，climatology 用验证集均值。
    """
    tot, n = 0.0, 0
    clim = None
    if baseline == "climatology":
        ys = [y for _, y, _ in loader]
        clim = torch.stack(ys).mean(dim=0, keepdim=True)
    for x, y, m in loader:
        pred = x if baseline == "persistence" else clim.expand_as(y)
        tot += masked_mse(pred, y, m).item() * x.shape[0]
        n += x.shape[0]
    return tot / max(n, 1)


@torch.no_grad()
def evaluate_model(model: nn.Module, loader: DataLoader) -> float:
    """Validation MSE of a trained model. 中文说明：已训练模型的验证集 MSE。"""
    model.eval()
    tot, n = 0.0, 0
    for x, y, m in loader:
        tot += masked_mse(model(x), y, m).item() * x.shape[0]
        n += x.shape[0]
    return tot / max(n, 1)


def main() -> None:
    ap = argparse.ArgumentParser(description="Train a tiny emulator on PlaSiC output")
    ap.add_argument("--data", required=True, help="dir with <var>.nc files")
    ap.add_argument("--vars", nargs="+", default=["tas", "ps", "pr"])
    ap.add_argument("--levels", type=int, default=3)
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--hidden", type=int, default=32)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--val-frac", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    full = PlasicDataset(args.data, args.vars, levels=args.levels)
    n_val = max(1, int(len(full) * args.val_frac))
    # chronological split: validate on the most recent months
    # 中文说明：按时间顺序划分，验证集取最近的月份（避免未来信息泄漏）
    n_train = len(full) - n_val
    train_ds = torch.utils.data.Subset(full, range(n_train))
    val_ds = torch.utils.data.Subset(full, range(n_train, len(full)))
    train_loader = DataLoader(train_ds, batch_size=args.batch, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch)

    print(f"train pairs: {len(train_ds)}, val pairs: {len(val_ds)}, "
          f"vars: {args.vars}")
    for name in ("persistence", "climatology"):
        print(f"baseline[{name}] val MSE: {evaluate(name, val_loader):.6f}")

    model = TinyCNN(len(args.vars), hidden=args.hidden)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    for epoch in range(1, args.epochs + 1):
        model.train()
        tot, n = 0.0, 0
        for x, y, m in train_loader:
            opt.zero_grad()
            loss = masked_mse(model(x), y, m)
            loss.backward()
            opt.step()
            tot += loss.item() * x.shape[0]
            n += x.shape[0]
        model.eval()
        val = evaluate_model(model, val_loader)
        print(f"epoch {epoch:3d} train MSE {tot / n:.6f} | val MSE {val:.6f}")

    out = Path(args.data) / "emulator.pt"
    torch.save({"state_dict": model.state_dict(), "vars": args.vars,
                "stats": full.stats}, out)
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
