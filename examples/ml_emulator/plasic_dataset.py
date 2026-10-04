#!/usr/bin/env python3
"""PyTorch Dataset for PlaSiC monthly-mean NetCDF output.

中文说明：读取 PlaSiC 月平均 NetCDF 输出的 PyTorch Dataset。
PlaSiC 把每个变量存为独立的 `<name>.nc` 文件（见 src/runtime/monthly_netcdf.c
的变量表，如 tas/ps/pr），维度为 (time, lat, lon) 或 (time, level, lat, lon)，
时间单位 "days since 1850-01-01"，掩膜点用 _FillValue 标记。

English: each variable lives in its own `<name>.nc` file with dimensions
(time, lat, lon) or (time, level, lat, lon). This dataset builds (X_t, X_{t+1})
next-step pairs with per-variable standardization; masked (_FillValue) points
are excluded from the loss via the returned mask.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import netCDF4
except ImportError as exc:  # pragma: no cover
    raise ImportError("netCDF4 is required: pip install netCDF4") from exc

try:
    import torch
    from torch.utils.data import Dataset
except ImportError as exc:  # pragma: no cover
    raise ImportError("torch is required: pip install torch") from exc


def load_variable(data_dir: str | Path, name: str,
                  levels: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Load one PlaSiC variable from ``<name>.nc``.

    中文说明：从 `<name>.nc` 读取单个变量，返回 (data, valid_mask)；
    data 形状为 (time, lat, lon)，3-D 变量可选取低 levels 层做垂直平均。

    Returns (data, valid_mask): data shaped (time, lat, lon) float32 with
    NaN at masked points; valid_mask bool array, True where data is valid.
    """
    path = Path(data_dir) / f"{name}.nc"
    if not path.exists():
        raise FileNotFoundError(f"PlaSiC output not found: {path}")
    with netCDF4.Dataset(path, "r") as ds:
        var = ds.variables[name]
        arr = var[:]  # netCDF4 auto-masks _FillValue -> MaskedArray
        if isinstance(arr, np.ma.MaskedArray):
            valid = ~np.ma.getmaskarray(arr)
            data = arr.filled(np.nan).astype(np.float64)
        else:
            data = np.asarray(arr, dtype=np.float64)
            fill = getattr(var, "_FillValue", None)
            if fill is not None and not np.isnan(fill):
                valid = data != fill
                data = np.where(valid, data, np.nan)
            else:
                valid = np.ones(data.shape, dtype=bool)
        if data.ndim == 4:  # (time, level, lat, lon)
            nlev = data.shape[1] if levels is None else min(levels, data.shape[1])
            # average valid levels only; a column is valid if any level is valid
            # 中文说明：多层取平均；只要有一层有效该列即有效
            with np.errstate(invalid="ignore"):
                data = np.nanmean(data[:, :nlev], axis=1)
            valid = valid[:, :nlev].any(axis=1)
        elif data.ndim != 3:
            raise ValueError(f"unexpected dims for {name}: {var.dimensions}")
    return data.astype(np.float32), valid


def compute_stats(data: np.ndarray) -> tuple[float, float]:
    """Mean/std over valid (non-NaN) points. 中文说明：只在有效点上算均值方差。"""
    return float(np.nanmean(data)), float(np.nanstd(data)) + 1e-8


class PlasicDataset(Dataset):
    """Next-step pairs (X_t, X_{t+1}) from PlaSiC monthly output.

    中文说明：把 PlaSiC 月平均输出组织成 next-step 样本对：
    输入 X_t 为 t 时刻的多变量场，目标为 t+1 时刻的同一组变量。
    每个变量独立标准化；返回的 mask 标出有效格点（loss 里只算这些点）。

    Each sample: x (C, H, W) standardized state at month t,
    y (C, H, W) standardized state at month t+1,
    m (C, H, W) bool mask of valid grid points.
    """

    def __init__(self, data_dir: str | Path, variables: list[str],
                 levels: int | None = 3,
                 stats: dict[str, tuple[float, float]] | None = None,
                 t_start: int = 0, t_end: int | None = None) -> None:
        self.variables = variables
        fields, masks = [], []
        for v in variables:
            data, valid = load_variable(data_dir, v, levels=levels)
            fields.append(data)
            # a grid point is valid only if valid at both t and t+1
            masks.append(valid)
        self.n_time = fields[0].shape[0]
        for f in fields:
            assert f.shape[0] == self.n_time, "time dimension mismatch"
        t_end = self.n_time - 1 if t_end is None else min(t_end, self.n_time - 1)
        self.t0, self.t1 = t_start, t_end
        self.stats = stats or {v: compute_stats(f) for v, f in zip(variables, fields)}
        self.fields = np.stack([
            (f - self.stats[v][0]) / self.stats[v][1]
            for v, f in zip(variables, fields)
        ])  # (C, T, H, W)
        # Zero-fill masked points; the mask (not the value) marks validity,
        # so NaN * 0 can never poison the masked loss.
        # 中文说明：掩膜点填 0，有效性由 mask 标识，避免 NaN 污染 masked loss。
        self.fields = np.nan_to_num(self.fields, nan=0.0).astype(np.float32)
        self.masks = np.stack(masks)  # (C, T, H, W)

    def __len__(self) -> int:
        return max(0, self.t1 - self.t0)

    def __getitem__(self, i: int):
        t = self.t0 + i
        x = torch.from_numpy(self.fields[:, t])
        y = torch.from_numpy(self.fields[:, t + 1])
        m = torch.from_numpy(self.masks[:, t + 1] & self.masks[:, t])
        return x, y, m
