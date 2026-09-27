#ifndef PLASIC_FLUXMOD_KERNELS_H
#define PLASIC_FLUXMOD_KERNELS_H

#include <stdint.h>

/*
 * C numerical kernels for fluxmod.f90.
 *
 * Every array is contiguous float storage.  Two-dimensional fields retain the
 * model's column-major layout, so element (horizontal, level) is stored at
 * horizontal + level * nhor.  All indices are zero based.
 */

void flux_surface_exchange_coefficients(
    int32_t nhor, int32_t low_level, int32_t surface_level, int32_t ntsa,
    float vonkarman, float kap, float sigma_low, float rdbrv, float gascon,
    float gravity, float zumin, float vdiff_b, float vdiff_c, float vdiff_d,
    const float *dt, const float *dq, const float *du, const float *dv,
    const float *dz0, const float *dls, float *dtsa, float *dtransh,
    float *dtransm, float *zabsu2, float *znl, float *zri, float *zrifh,
    float *zrifm);

void flux_mkstress(
    int32_t nhor, int32_t low_level, int32_t surface_level, int32_t ndheat,
    float gravity, float deltsec2, float gascon, float dsigma_low, float cpd,
    float cpv_cpd_minus1, const float *dtransm, const float *dt, const float *dq,
    const float *du, const float *dv, const float *dp, float *dudt,
    float *dvdt, float *dtdt, float *dtaux, float *dtauy, float *dust3,
    float *zun, float *zvn);

void flux_mkshfl(
    int32_t nhor, int32_t low_level, int32_t surface_level, float kap,
    float sigma_low, float gravity, float deltsec2, float gascon,
    float dsigma_low, float cpd, float cpv_cpd_minus1, const float *dtransh,
    const float *dt, const float *dq, const float *dp, float *dtdt,
    float *dshfl, float *dshdt, float *dtsa, float *ztn, float *zfac);

void flux_mkevap(
    int32_t nhor, int32_t low_level, int32_t surface_level, float gravity,
    float deltsec, float deltsec2, float gascon, float dsigma_low, float tmelt,
    float lv, float ls, float ra2, float ra4, const float *dtransh,
    const float *dt, const float *dq, const float *dp, const float *dwetfac,
    const float *dls, const float *dwatc, float *dqdt, float *devap,
    float *dlhfl, float *dlhdt, float *zqn);

#endif
