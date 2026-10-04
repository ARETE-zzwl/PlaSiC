#!/usr/bin/env python3
"""Train a tiny next-step emulator on PlaSiC output.

在 PlaSiC 输出上训练小型 next-step emulator 的示例：持续性/气候态
两个基线加一个小 CNN，只在掩膜外的有效格点上算 MSE。先拿合成数据
跑通，再换成真实的 PlaSiC 输出。

Baselines (persistence, climatology) plus a small CNN emulator,
chronological train/val split, masked MSE. Run on synthetic data
first, then on real PlaSiC output.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset

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
        # residual: predict the increment, not the full state
        return x + self.net(x)


def masked_mse(pred: torch.Tensor, target: torch.Tensor,
               mask: torch.Tensor) -> torch.Tensor:
    se = (pred - target) ** 2
    return (se * mask).sum() / mask.sum().clamp(min=1)


@torch.no_grad()
def evaluate(predictor, loader: DataLoader,
             device: torch.device | None = None) -> float:
    """Validation MSE of a baseline or a model.

    predictor 为 'persistence' / 'climatology' / nn.Module。
    climatology 取验证集的时间平均。
    """
    if predictor == "climatology":
        ys = torch.cat([y for _, y, _ in loader], dim=0)
        clim = ys.mean(dim=0, keepdim=True)  # (1, C, H, W)
    tot, n = 0.0, 0
    for x, y, m in loader:
        if device is not None:
            x, y, m = x.to(device), y.to(device), m.to(device)
        if predictor == "persistence":
            pred = x
        elif predictor == "climatology":
            pred = clim.to(x.device).expand_as(y)
        else:
            predictor.eval()
            pred = predictor(x)
        tot += masked_mse(pred, y, m).item() * x.shape[0]
        n += x.shape[0]
    return tot / max(n, 1)


@torch.no_grad()
def per_variable_mse(model: nn.Module, loader: DataLoader,
                     variables: list[str],
                     device: torch.device | None = None) -> dict[str, float]:
    """One-step MSE per variable. 分变量的单步 MSE。"""
    model.eval()
    tot = np.zeros(len(variables))
    cnt = np.zeros(len(variables))
    for x, y, m in loader:
        if device is not None:
            x, y, m = x.to(device), y.to(device), m.to(device)
        se = (model(x) - y) ** 2  # (B, C, H, W)
        for c in range(len(variables)):
            mc = m[:, c]
            tot[c] += (se[:, c] * mc).sum().item()
            cnt[c] += mc.sum().item()
    return {v: tot[c] / max(cnt[c], 1) for c, v in enumerate(variables)}


@torch.no_grad()
def rollout_mse(model: nn.Module, dataset: PlasicDataset,
                i_start: int, i_end: int, steps: int = 4,
                device: torch.device | None = None) -> list[float]:
    """Autoregressive rollout MSE, one value per lead step.

    自回归多步 rollout：把模型输出喂回去当下一步输入，
    看误差随 lead time 怎么涨。这是 emulator 的常规评估。
    """
    model.eval()
    errs = np.zeros(steps)
    n = 0
    for i in range(i_start, min(i_end, len(dataset) - steps)):
        x, _, _ = dataset[i]
        cur = x.unsqueeze(0)
        if device is not None:
            cur = cur.to(device)
        for s in range(1, steps + 1):
            cur = model(cur)
            _, y_true, m_true = dataset[i + s - 1]
            if device is not None:
                y_true, m_true = y_true.to(device), m_true.to(device)
            errs[s - 1] += masked_mse(cur, y_true.unsqueeze(0),
                                      m_true.unsqueeze(0)).item()
        n += 1
    return (errs / max(n, 1)).tolist()


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
    ap.add_argument("--device", default="auto",
                    help="'auto', 'cpu', or 'cuda'")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device(
        "cuda" if (args.device == "cuda" or
                   (args.device == "auto" and torch.cuda.is_available()))
        else "cpu")
    print(f"device: {device}")

    full = PlasicDataset(args.data, args.vars, levels=args.levels)
    n_val = max(1, int(len(full) * args.val_frac))
    n_train = len(full) - n_val
    # Split chronologically: validate on the most recent months, so no
    # future information leaks into training.
    train_loader = DataLoader(Subset(full, range(n_train)),
                              batch_size=args.batch, shuffle=True)
    val_loader = DataLoader(Subset(full, range(n_train, len(full))),
                            batch_size=args.batch)

    print(f"train pairs: {n_train}, val pairs: {n_val}, vars: {args.vars}")
    print(f"baseline[persistence] val MSE: {evaluate('persistence', val_loader):.6f}")
    print(f"baseline[climatology] val MSE: {evaluate('climatology', val_loader):.6f}")

    model = TinyCNN(len(args.vars), hidden=args.hidden).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    history = []
    for epoch in range(1, args.epochs + 1):
        model.train()
        tot, n = 0.0, 0
        for x, y, m in train_loader:
            x, y, m = x.to(device), y.to(device), m.to(device)
            opt.zero_grad()
            loss = masked_mse(model(x), y, m)
            loss.backward()
            opt.step()
            tot += loss.item() * x.shape[0]
            n += x.shape[0]
        val = evaluate(model, val_loader, device)
        history.append({"epoch": epoch, "train": tot / n, "val": val})
        print(f"epoch {epoch:3d}  train {tot / n:.6f}  val {val:.6f}")

    out = Path(args.data) / "emulator.pt"
    torch.save({"state_dict": model.state_dict(), "vars": args.vars,
                "stats": full.stats}, out)
    print(f"saved -> {out}")

    # Final report: per-variable one-step MSE + rollout error growth.
    # 最终报告：分变量单步 MSE，以及 rollout 误差随 lead time 的增长。
    per_var = per_variable_mse(model, val_loader, args.vars, device)
    steps = 4
    ro = rollout_mse(model, full, n_train, len(full), steps=steps, device=device)
    print("per-variable val MSE:",
          {v: f"{e:.6f}" for v, e in per_var.items()})
    print("rollout val MSE by lead month:",
          {f"+{s}m": f"{e:.6f}" for s, e in enumerate(ro, 1)})

    metrics = {
        "vars": args.vars,
        "n_train": n_train,
        "n_val": n_val,
        "history": history,
        "per_variable_val_mse": per_var,
        "rollout_val_mse": {f"+{s}m": e for s, e in enumerate(ro, 1)},
    }
    mpath = Path(args.data) / "metrics.json"
    mpath.write_text(json.dumps(metrics, indent=2))
    print(f"metrics -> {mpath}")


if __name__ == "__main__":
    main()
