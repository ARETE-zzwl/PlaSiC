# Offline ERA5 -> PlaSiC spectral helper.
#
# The helper must be compiled with the same atmospheric grid and level count as
# the model executable that consumes its restart files.  Defaults match the
# forecast workflow (T85/L25 with the slab ocean); every setting can be
# overridden on the make command line.  For example:
#
#   make -C experiments/forecast_era5/tools \
#        NLAT=128 NLEV=40 OCEAN_MODEL=3d DEEP_OCEAN_LEVELS=12 \
#        MODEL_BUILD=../../../build/T85_L40_O12
#
# 三维海洋不是辅助工具的必选项；只有显式设置 ``OCEAN_MODEL=3d`` 才会打开，
# ``DEEP_OCEAN_LEVELS`` 也可以按实际模式构建自由设置，不再固定为 16 层。

PROJECT_ROOT ?= ../../..
SHTNS_PREFIX ?= $(HOME)/.local/shtns-3.7.5

# Keep these defaults in sync with common.py.  NLAT determines the triangular
# truncation exactly as in src/Makefile; NLEV and ocean settings are independent.
NLAT ?= 128
NLEV ?= 25
OCEAN_MODEL ?= slab
DEEP_OCEAN_LEVELS ?= 16
TRUNCATION ?= $(shell expr \( $(NLAT) \* 2 - 1 \) / 3)

CC ?= cc
CFLAGS ?= -O2
CPPFLAGS += -std=c11 -I$(PROJECT_ROOT)/src -I$(PROJECT_ROOT)/src/state \
	-DPLASIC_NLAT=$(NLAT) -DPLASIC_NLEV=$(NLEV) -DPLASIC_NPRO=1 \
	-DPLASIC_USE_MPI=0 -DPLASIC_DEEP_NLEV=$(DEEP_OCEAN_LEVELS) \
	-DPLASIC_DEEP_OCEAN=$(if $(filter 3d,$(OCEAN_MODEL)),1,0)

# Override MODEL_BUILD when the model was built into a custom directory.  The
# default follows src/Makefile's serial-build naming convention.
OCEAN_BUILD_SUFFIX = $(if $(filter 3d,$(OCEAN_MODEL)),_O$(DEEP_OCEAN_LEVELS),)
MODEL_BUILD ?= $(PROJECT_ROOT)/build/T$(TRUNCATION)_L$(NLEV)$(OCEAN_BUILD_SUFFIX)
LDLIBS = $(MODEL_BUILD)/libplasic.a -L$(SHTNS_PREFIX)/lib -lshtns \
	-L$(shell brew --prefix fftw)/lib -lfftw3 \
	-L$(shell brew --prefix netcdf)/lib -lnetcdf \
	-framework Accelerate -lm

all: plasic_era5_tools

plasic_era5_tools: plasic_era5_tools.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f plasic_era5_tools

.PHONY: all clean
