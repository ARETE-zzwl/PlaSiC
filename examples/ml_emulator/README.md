# ML Emulator Example — train a neural network on PlaSiC output

[中文说明](#中文说明) · [English](#english)

---

## English

### What is this?

A minimal, self-contained example showing **machine-learning researchers**
(audience #2 of PlaSiC) how to train a neural-network **emulator** on PlaSiC
model output. The emulator learns next-step prediction:

```
X(t) = [tas, ps, pr, ...](t)   →   X̂(t+Δt)
```

i.e. given the model state at month *t*, predict the state at month *t+1*.
This is the standard "learn the dynamics" setup used in ML-for-weather
research (cf. NeuralGCM, GraphCast-style next-step training).

### Files

| File | Purpose |
|---|---|
| `plasic_dataset.py` | PyTorch `Dataset` reading PlaSiC monthly NetCDF output (`<var>.nc` files), building `(X_t, X_{t+1})` pairs with per-variable normalization |
| `train_emulator.py` | Training script: persistence & climatology baselines + a small CNN emulator, train/val split, MSE reporting |
| `make_synthetic_sample.py` | Generates tiny synthetic `<var>.nc` files mimicking PlaSiC output so the example runs without a full model run |
| `requirements.txt` | Python dependencies |

### Quickstart (no model run needed)

```bash
pip install -r requirements.txt

# 1. create synthetic PlaSiC-like output (tas/ps/pr, 24 months, T21-ish grid)
python make_synthetic_sample.py --out data/synthetic

# 2. train the emulator on it
python train_emulator.py --data data/synthetic --epochs 10
```

Expected: the CNN beats the persistence baseline on validation MSE, because
the synthetic data contains a predictable propagating signal.

### Using real PlaSiC output

1. Run any PlaSiC experiment with monthly NetCDF output enabled (e.g. the
   historical experiment, tutorial §6.2). Each variable is written to
   `<name>.nc` (`tas.nc`, `ps.nc`, `pr.nc`, …) with CF conventions —
   see `src/runtime/monthly_netcdf.c` for the full variable table.
2. Point the scripts at the output directory:
   ```bash
   python train_emulator.py --data /path/to/plasic/output --vars tas ps pr \
       --epochs 50 --hidden 64
   ```
3. For 3-D variables (`ta`, `ua`, `va`, `hus`), select levels with
   `--levels 3` (uses the lowest 3 levels) or extend `plasic_dataset.py`.

### Notes / limitations

- Monthly-mean data is coarse for weather emulation; this example is a
  **template**, not a production emulator. For serious work, use higher-
  frequency output or the restart files (tutorial §6.4).
- `_FillValue` (land/ocean-masked points) is masked out of the loss.
- Contributions welcome: additional architectures, proper train/val/test
  splits by year, rollout (multi-step) evaluation.

---

## 中文说明

### 这是什么？

一个最小可运行的示例，展示**机器学习研究者**（PlaSiC 的第二类目标用户）
如何在 PlaSiC 模式输出上训练神经网络** emulator**。Emulator 学习单步预测：

```
X(t) = [tas, ps, pr, ...](t)   →   X̂(t+Δt)
```

即给定第 *t* 个月的模式状态，预测第 *t+1* 个月。这是 ML-for-weather
研究的标准 "learn the dynamics" 范式（参见 NeuralGCM、GraphCast 式的
next-step 训练）。

### 文件说明

| 文件 | 用途 |
|---|---|
| `plasic_dataset.py` | PyTorch `Dataset`，读取 PlaSiC 月平均 NetCDF 输出（`<var>.nc` 文件），构造 `(X_t, X_{t+1})` 样本对，逐变量标准化 |
| `train_emulator.py` | 训练脚本：持续性预报 & 气候态基线 + 小型 CNN emulator，训练/验证划分，MSE 评估 |
| `make_synthetic_sample.py` | 生成模仿 PlaSiC 输出的小型合成 `<var>.nc` 文件，无需真实跑模式即可运行示例 |
| `requirements.txt` | Python 依赖 |

### 快速开始（无需跑模式）

```bash
pip install -r requirements.txt

# 1. 生成合成的 PlaSiC 风格输出（tas/ps/pr，24 个月，T21 量级网格）
python make_synthetic_sample.py --out data/synthetic

# 2. 在上面训练 emulator
python train_emulator.py --data data/synthetic --epochs 10
```

预期结果：CNN 在验证集 MSE 上打败持续性基线，因为合成数据里含有一个
可预测的传播信号。

### 使用真实的 PlaSiC 输出

1. 运行任意开启月平均 NetCDF 输出的 PlaSiC 试验（如 historical 试验，
   教程 §6.2）。每个变量存为 `<name>.nc`（`tas.nc`、`ps.nc`、`pr.nc`…），
   遵循 CF 规范——完整变量表见 `src/runtime/monthly_netcdf.c`。
2. 把脚本指向输出目录：
   ```bash
   python train_emulator.py --data /path/to/plasic/output --vars tas ps pr \
       --epochs 50 --hidden 64
   ```
3. 三维变量（`ta`、`ua`、`va`、`hus`）可用 `--levels 3` 取低 3 层，
   或自行扩展 `plasic_dataset.py`。

### 说明与局限

- 月平均资料对天气 emulator 来说时间分辨率较粗；本示例是**模板**，
  不是生产级 emulator。严肃研究请用更高频输出或重启文件（教程 §6.4）。
- `_FillValue`（陆/海掩膜点）在 loss 中被 mask 掉。
- 欢迎贡献：更多网络结构、按年划分的 train/val/test、多步 rollout 评估。
