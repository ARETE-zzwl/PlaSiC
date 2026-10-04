#!/usr/bin/env python3
"""Train a tiny next-step emulator on PlaSiC output.

中文说明：在 PlaSiC 输出上训练小型 next-step emulator 的示例。
含持续性/气候态两个基线和一个小 CNN；只在掩膜外的有效格点上算 MSE。
先用合成数据跑通，再换成真实的 PlaSiC 输出。

Baselines (persistence, climatology) + a small CNN emulator,
chronological train/val split, masked MSE. Run on synthetic data
first, then on real PlaSiC output.
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
    def __init__(self, channels: int, hidden: int = 32) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(channels, hidden, 3, padding=1), nn.ReLU(),
            nn.Conv2d(hidden, hidden, 3, padding=1), nn.ReLU(),
            nn.Conv2d(hidden, channels, 3, padding=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)  # residual: predict the increment, not the state


def masked_mse(pred: torch.Tensor, target: torch.Tensor,
               mask: torch.Tensor) -> torch.Tensor:
    se = (pred - target) ** 2
    return (se * mask).sum() / mask.sum().clamp(min=1)


@torch.no_grad()
def evaluate(predictor, loader: DataLoader) -> float:
    """MSE of a predictor on a loader.

    predictor: 'persistence' | 'climatology' | nn.Module.
    中文说明：基线或模型的验证 MSE；climatology 取验证集时间平均。
    """
    if predictor == "climatology":
        ys = torch.cat([y for _, y, _ in loader], dim=0)
        clim = ys.mean(dim=0, keepdim=True)  # (1, C, H, W)
    tot, n = 0.0, 0
    for x, y, m in loader:
        if predictor == "persistence":
            pred = x
        elif predictor == "climatology":
            pred = clim.expand_as(y)
        else:
            predictor.eval()
            pred = predictor(x)
        tot += masked_mse(pred, y, m).item() * x.shape[0]
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
    n_train = len(full) - n_val
    # Chronological split: validate on the most recent months so no
    # future information leaks into training.
    train_loader = DataLoader(torch.utils.data.Subset(full, range(n_train)),
                              batch_size=args.batch, shuffle=True)
    val_loader = DataLoader(torch.utils.data.Subset(full, range(n_train, len(full))),
                            batch_size=args.batch)

    print(f"train pairs: {n_train}, val pairs: {n_val}, vars: {args.vars}")
    print(f"baseline[persistence] val MSE: {evaluate('persistence', val_loader):.6f}")
    print(f"baseline[climatology] val MSE: {evaluate('climatology', val_loader):.6f}")

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
        val = evaluate(model, val_loader)
        print(f"epoch {epoch:3d}  train {tot / n:.6f}  val {val:.6f}")

    out = Path(args.data) / "emulator.pt"
    torch.save({"state_dict": model.state_dict(), "vars": args.vars,
                "stats": full.stats}, out)
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
