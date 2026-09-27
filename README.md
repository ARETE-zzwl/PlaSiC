# PlaSiC

PlaSiC（Planet Simulator in C）是一套中等复杂度的大气环流模式。PlaSiC 完全用 C11 重写，不再依赖 Fortran。

## 构建依赖

PlaSiC 的水平标量和矢量球谐变换由 [SHTns](https://nschaeff.bitbucket.io/shtns/) 完成，SHTns 内部使用 FFTW3。本项目已使用 SHTns 3.7.5 测试。

下面是在 macOS 上将 SHTns 安装到用户持久目录的完整流程。安装位置为 `$HOME/.local/shtns-3.7.5`。

### 1. 安装基础依赖

首先安装 [Homebrew](https://brew.sh/)，然后安装 FFTW3、NetCDF-C 和
`pkg-config`：

```bash
brew install fftw netcdf pkg-config
```

如果需要编译和运行 MPI 版本，还需要 Open MPI：

```bash
brew install open-mpi
```

### 2. 设置版本、安装目录和编译目录

```bash
SHTNS_VERSION="3.7.5"
SHTNS_PREFIX="$HOME/.local/shtns-${SHTNS_VERSION}"
SHTNS_WORK="$HOME/Downloads/shtns-${SHTNS_VERSION}-build"

mkdir -p "$SHTNS_WORK"
cd "$SHTNS_WORK"
```

`SHTNS_PREFIX` 是最终安装目录；`SHTNS_WORK` 只用于保存下载的源码和中间编译文件。

### 3. 下载并校验 SHTns 3.7.5

```bash
curl -fL \
  "https://files.pythonhosted.org/packages/35/60/241291c0a737d269efd0ae02a488172721f1ad4b226729613cf8d12d3a4f/shtns-3.7.5.tar.gz" \
  -o "shtns-3.7.5.tar.gz"
```

### 4. 解压并配置

```bash
tar -xzf "shtns-3.7.5.tar.gz"
cd "shtns-3.7.5"

FFTW_PREFIX="$(brew --prefix fftw)"

CPPFLAGS="-I${FFTW_PREFIX}/include" \
LDFLAGS="-L${FFTW_PREFIX}/lib" \
./configure \
  --prefix="$SHTNS_PREFIX" \
  --disable-openmp
```

当前 PlaSiC 的 MPI 教学路径由每个进程各自执行完整 SHTns 变换。因此这里默认关闭 SHTns 内部 OpenMP，避免 MPI 进程和 OpenMP 线程同时抢占 CPU。

### 5. 编译 SHTns

```bash
make -j"$(sysctl -n hw.logicalcpu)"
```

编译完成后应该生成静态库：

```bash
ls -lh libshtns.a
```

### 6. 安装到持久目录

SHTns 3.7.5 的源码包在非 GPU 构建下执行 `make install` 时，可能仍会尝试复制源码包中不存在的 CUDA 头文件。PlaSiC 只需要 `shtns.h` 和 `libshtns.a`，因此直接安装这两个文件：

```bash
mkdir -p "$SHTNS_PREFIX/include"
mkdir -p "$SHTNS_PREFIX/lib"

install -m 0644 shtns.h \
  "$SHTNS_PREFIX/include/shtns.h"

install -m 0644 libshtns.a \
  "$SHTNS_PREFIX/lib/libshtns.a"
```

检查安装结果：

```bash
ls -lh \
  "$SHTNS_PREFIX/include/shtns.h" \
  "$SHTNS_PREFIX/lib/libshtns.a"
```

### 7. 让 PlaSiC 长期找到 SHTns

在 `~/.zshrc` 中加入：

```bash
export SHTNS_PREFIX="$HOME/.local/shtns-3.7.5"
```

然后重新加载 shell 配置：

```bash
source ~/.zshrc
```

确认环境变量：

```bash
echo "$SHTNS_PREFIX"
```

应该输出类似：

```text
/Users/your-name/.local/shtns-3.7.5
```

PlaSiC 的 Makefile 会读取这个环境变量，以后构建时不需要每次手动传入 `SHTNS_PREFIX=...`。如果不想修改 `~/.zshrc`，也可在每次构建时显式指定：

```bash
make -C src SHTNS_PREFIX="$HOME/.local/shtns-3.7.5"
```

### 8. 编译和测试 PlaSiC

进入 PlaSiC 项目根目录，先运行严格警告构建：

```bash
make -C src strict
```

编译串行版 PlaSiC：

```bash
make -C src
```

编译 MPI 双进程版 PlaSiC：

```bash
make -C src MPI=1 NPRO=2
```

Makefile 会优先使用 `pkg-config` 定位 FFTW3，也可直接覆盖 `SHTNS_CFLAGS`/`SHTNS_LIBS` 和 `FFTW_CFLAGS`/`FFTW_LIBS`。当前安装的是静态 `libshtns.a`，因此运行 PlaSiC 时不需要另外设置 `DYLD_LIBRARY_PATH`。


## 在线文档

文档支持本地预览（需要先创建虚拟环境）：

```bash
cd docs-site
python3 -m venv .venv
.venv/bin/pip install mkdocs-material
.venv/bin/mkdocs serve -f mkdocs.yml
```

打开 <http://127.0.0.1:8000/PLASIM/> 即可浏览。

文档站目前支持在线文档，登录 https://sunmoumou1.github.io/PlaSiC/tutorial/

## 开源许可与上游署名

Except for third-party materials stated otherwise, PlaSiC's C code, build and test scripts, desktop application, this online tutorial, and the website source are released under the **GNU General Public License version 3 or later** (SPDX: `GPL-3.0-or-later`).

The GPL permits research, teaching, modification, redistribution, and commercial use. When distributing binaries or derivative works governed by the GPL, distributors must provide the complete corresponding source code, retain applicable notices, and release GPL-covered derivative works under compatible GPL terms. This software and documentation are provided without any warranty.

Third-party materials and external dependencies remain subject to their own licenses and are not relicensed by being included in, or used with, this project.

PlaSiC is a derivative of [PlaSim (Planet Simulator)](https://github.com/HartmutBorth/PLASIM).For scientific attribution, the upstream model description should be cited as Fraedrich et al. (2005; full reference below). Upstream PlaSim is licensed under `GPL-2.0-or-later`, and the original author's copyright and license notices remain in force; a copy of the GPLv2 text is included. Since 2026, PlaSiC has reimplemented the model in C11, restructured the software, removed the SimBA module, do the historical & 4xCO2 & 1pctCO2 experiments and added the desktop application and this online tutorial. The combined PlaSiC work is released under `GPL-3.0-or-later`, exercising the "GPL v2 or any later version" permission granted by the upstream project.

This project owes a great deal to the teaching materials and open-source work of others. In particular, *Introduction to Climate Modelling* by Stocker (2011) and the [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) project (Klöwer et al., 2024) guided much of the design of PlaSiC. The former shaped how the physics and mathematics are presented, while the latter is an inspiring model of an interactive, approachable atmospheric GCM. Without these two resources, PlaSiC could not have been designed. The original authors are gratefully acknowledged.

- Stocker, T. (2011). *Introduction to Climate Modelling*. Springer Science & Business Media.
- Klöwer, M., Gelbrecht, M., Hotta, D., Willmert, J., Silvestri, S., Wagner, G. L., White, A., Hatfield, S., Kimpson, T., Constantinou, N. C., & Hill, C. (2024). SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity and extensibility. *Journal of Open Source Software*, *9*(98), 6323. [https://doi.org/10.21105/joss.06323](https://doi.org/10.21105/joss.06323)


