# Offline ERA5 -> PlaSiC T85/L25 helper.
#
# Links against the model's own kernels so the spherical-harmonic analysis is
# identical to the model's.  Build serial T85/L25 objects first:
#   make -C src all MPI=0 NPRO=1 NLAT=128 NLEV=25

PROJECT_ROOT ?= ../../..
SHTNS_PREFIX ?= $(HOME)/.local/shtns-3.7.5

CC ?= cc
CFLAGS ?= -O2
CPPFLAGS += -std=c11 -I$(PROJECT_ROOT)/src -I$(PROJECT_ROOT)/src/state \
	-DPLASIC_NLAT=128 -DPLASIC_NLEV=25 -DPLASIC_NPRO=1 \
	-DPLASIC_USE_MPI=0 -DPLASIC_DEEP_OCEAN=1 -DPLASIC_DEEP_NLEV=16

MODEL_BUILD = $(PROJECT_ROOT)/build/T85_L25_O16
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
