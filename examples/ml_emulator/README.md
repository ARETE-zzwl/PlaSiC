# ML Emulator Example — train a neural network on PlaSiC output

[中文说明](#中文说明) · [English](#english)

---

## English

### What is this?

A small working example: train a neural net to do next-step prediction on
PlaSiC output. Given the model state at month *t*, predict month *t+1*:

```
X(t) = [tas, ps, pr, ...](t)   →   X̂(t+Δt)
```

Same idea as NeuralGCM / GraphCast-style next-step training, just much
smaller. If you come from ML and want to poke at PlaSiC data, start here.

### Files

| File | Purpose |
|---|---|
| `plasic_dataset.py` | PyTorch `Dataset` reading PlaSiC monthly NetCDF output (`<var>.nc` files), building `(X_t, X_{t+1})` pairs with per-variable normalization |
| `train_emulator.py` | Training script: persistence & climatology baselines + a small CNN emulator, chronological train/val split, masked MSE, per-variable metrics, multi-step rollout evaluation; writes `emulator.pt` and `metrics.json` |
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

Useful flags: `--device cuda` (auto-detects by default), `--hidden 64`,
`--epochs 50`, `--lr 3e-4`.

### Reading the output

A typical run prints something like:

```
device: cpu
train pairs: 17, val pairs: 6, vars: ['tas', 'ps', 'pr']
baseline[persistence] val MSE: 0.512340
baseline[climatology] val MSE: 0.498112
epoch   1  train 0.501233  val 0.495012
...
epoch  10  train 0.201455  val 0.223311
saved -> data/synthetic/emulator.pt
per-variable val MSE: {'tas': '0.039...', 'ps': '0.004...', 'pr': '0.650...'}
rollout val MSE by lead month: {'+1m': '0.24...', '+2m': '0.31...', '+3m': '0.39...', '+4m': '0.47...'}
metrics -> data/synthetic/metrics.json
```

What to look at:

- **CNN vs baselines.** If the CNN can't beat persistence on one-step MSE,
  something is off — check the data first, not the model.
- **Per-variable MSE.** `pr` (precipitation) is usually the worst by far.
  That's normal: it's intermittent and skewed, not a bug in your setup.
- **Rollout error growth.** Error should grow with lead time, roughly
  monotonically. If `+4m` is *better* than `+1m`, your rollout is probably
  leaking (e.g. teacher forcing sneaked in).
- **`metrics.json`** holds the per-epoch history plus the final report, so
  you can plot learning curves without re-running.

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
3. For 3-D variables (`ta`, `ua`, `va`, `hus`), `--levels 3` averages the
   lowest 3 levels; edit `plasic_dataset.py` if you want something else
   (e.g. a single level, or all levels as extra channels).

### Try it yourself

- Swap `TinyCNN` for a two-layer MLP and compare: where does the CNN win,
  and why?
- Train on anomalies instead of full fields (subtract the climatology in
  `plasic_dataset.py`). Does the rollout stay stable longer?
- Push `--epochs` until validation MSE stops improving, then check whether
  the rollout at `+4m` also stopped improving. They don't always agree.
- Add `pr` later than `tas`/`ps`: how much does precipitation hurt or help
  the other variables?

### FAQ

**NaN loss in the first epoch?**
Almost always a normalization problem: a variable with zero variance
(e.g. a constant field) gives `std = 0`. `compute_stats` guards with
`+ 1e-8`, but check your input files if it still happens.

**`FileNotFoundError: PlaSiC output not found`?**
The script expects one `.nc` file per variable, named exactly `<var>.nc`
inside `--data`. PlaSiC writes these names by default; if you renamed
them, rename back or adjust `load_variable`.

**Can I use this on the restart files instead of monthly means?**
Not directly — restarts are spectral coefficients, a different layout.
Monthly means are the intended input here; see tutorial §6.4 for what's
in the restart files.

**Why is rollout error so much larger than one-step error?**
Because rollout feeds the model's own (imperfect) predictions back in.
One-step error is optimistic; rollout is the honest number. This gap is
the main reason emulators are evaluated autoregressively.

### Notes / limitations

- Monthly means are coarse for weather emulation. This is a starting
  template — for real work you'd want higher-frequency output or the
  restart files (tutorial §6.4).
- `_FillValue` (land/ocean-masked points) is masked out of the loss.
- Things worth trying next: more architectures, train/val/test splits
  by year, proper early stopping.

---

## 中文说明

### 这是什么？

一个能跑通的小例子：在 PlaSiC 模式输出上训练神经网络做单步预测。
给定第 *t* 个月的模式状态，预测第 *t+1* 个月：

```
X(t) = [tas, ps, pr, ...](t)   →   X̂(t+Δt)
```

思路跟 NeuralGCM / GraphCast 那种 next-step 训练一样，只是小得多。
从 ML 过来的、想摆弄一下 PlaSiC 数据的人，可以从这儿入手。

### 文件说明

| 文件 | 用途 |
|---|---|
| `plasic_dataset.py` | PyTorch `Dataset`，读取 PlaSiC 月平均 NetCDF 输出（`<var>.nc` 文件），构造 `(X_t, X_{t+1})` 样本对，逐变量标准化 |
| `train_emulator.py` | 训练脚本：持续性预报 & 气候态基线 + 小型 CNN emulator，按时间划分训练/验证，mask MSE，分变量指标和多步 rollout 评估；输出 `emulator.pt` 和 `metrics.json` |
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

常用参数：`--device cuda`（默认自动检测）、`--hidden 64`、
`--epochs 50`、`--lr 3e-4`。

### 读懂输出

一次典型的运行会打印：

```
device: cpu
train pairs: 17, val pairs: 6, vars: ['tas', 'ps', 'pr']
baseline[persistence] val MSE: 0.512340
baseline[climatology] val MSE: 0.498112
epoch   1  train 0.501233  val 0.495012
...
epoch  10  train 0.201455  val 0.223311
saved -> data/synthetic/emulator.pt
per-variable val MSE: {'tas': '0.039...', 'ps': '0.004...', 'pr': '0.650...'}
rollout val MSE by lead month: {'+1m': '0.24...', '+2m': '0.31...', '+3m': '0.39...', '+4m': '0.47...'}
metrics -> data/synthetic/metrics.json
```

看这几个地方：

- **CNN 有没有打赢基线。** 单步 MSE 上如果连持续性预报都打不过，
  先查数据，别先调模型。
- **分变量 MSE。** `pr`（降水）通常差得最多，这正常：降水间歇性强、
  分布偏，不是你的配置有问题。
- **rollout 误差增长。** 误差应该随 lead time 单调上涨。如果 `+4m`
  反而比 `+1m` 好，多半是 rollout 写漏了（比如混进了 teacher forcing）。
- **`metrics.json`** 存了每轮的 train/val 和最终报告，画学习曲线
  不用重跑。

### 使用真实的 PlaSiC 输出

1. 运行任意开启月平均 NetCDF 输出的 PlaSiC 试验（如 historical 试验，
   教程 §6.2）。每个变量存为 `<name>.nc`（`tas.nc`、`ps.nc`、`pr.nc`…），
   遵循 CF 规范——完整变量表见 `src/runtime/monthly_netcdf.c`。
2. 把脚本指向输出目录：
   ```bash
   python train_emulator.py --data /path/to/plasic/output --vars tas ps pr \
       --epochs 50 --hidden 64
   ```
3. 三维变量（`ta`、`ua`、`va`、`hus`）用 `--levels 3` 取低 3 层做平均；
   想要别的处理（单层、或全层当通道）就改 `plasic_dataset.py`。

### 自己动手试试

- 把 `TinyCNN` 换成两层 MLP 对比：CNN 赢在哪，为什么？
- 改成预测距平而不是全场（在 `plasic_dataset.py` 里减去气候态），
  rollout 会不会更稳定？
- 把 `--epochs` 加到验证集 MSE 不再下降，看 `+4m` 的 rollout 是不是
  也不再变好——这两者不总是一致的。
- 先只用 `tas`/`ps` 训练，再加入 `pr`：降水对别的变量是帮忙还是添乱？

### 常见问题

**第一个 epoch 就 NaN？**
多半是标准化的问题：方差为 0 的变量（比如常数场）会导致除零。
`compute_stats` 里加了 `+ 1e-8` 兜底，还出现的话检查输入文件。

**`FileNotFoundError: PlaSiC output not found`？**
脚本要求 `--data` 里每个变量恰好对应一个 `<var>.nc` 文件。
PlaSiC 默认就是这么写的；改过名的改回去，或改 `load_variable`。

**能直接用重启文件（restart）代替月平均吗？**
不行，重启文件是谱系数，格式完全不同。这里认准月平均输出；
重启文件里有什么见教程 §6.4。

**为什么 rollout 误差比单步大这么多？**
因为 rollout 把模型自己的（不完美的）预测喂回去当输入。
单步误差是乐观估计，rollout 才是实在数。做 emulator 评估看
rollout，主要就是因为这个 gap。

### 说明与局限

- 月平均资料做天气 emulator 时间分辨率偏粗。这是个起点模板，
  真要做研究得用更高频的输出或重启文件（教程 §6.4）。
- `_FillValue`（陆/海掩膜点）在 loss 里被 mask 掉了。
- 下一步可以试试：更多网络结构、按年划分 train/val/test、
  正经的 early stopping。
