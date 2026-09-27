#ifndef PLASIC_SLAB_OCEAN_H
#define PLASIC_SLAB_OCEAN_H

#include <stdint.h>

/*
 * Numerical kernels for oceanmod.f90. Arrays use the model's column-major
 * storage. Month indices passed to the interpolation kernel are the original
 * 0..13 indices used by OCEANMOD.
 *
 * The mixed-layer ocean is forced to a single layer and has no vertical or
 * horizontal diffusion: each ocean point only exchanges heat with the
 * atmosphere, the sea ice and the deep ocean.
 */

void ocean_initialize_parameters(
    int32_t nhor, int32_t ntspd, float solar_day,
    float *taunc, float layer_depth, float *ymld, float *dtmix);

void ocean_interpolate_cycle(
    int32_t nhor, int32_t month1, int32_t month2, float weight,
    const float *monthly, float *field);

void ocean_copy_outputs(
    int32_t nhor, const float *ysst, const float *ymld,
    const float *yiflux, float *psst, float *pmld, float *piflux);

void ocean_copy_inputs(
    int32_t nhor, float crhoi, const float *picec, const float *piced,
    const float *pheat, const float *ppme, const float *ptaux,
    const float *ptauy, const float *pust3,
    const float *psnow, const float *pcliced, float *yicec, float *yiced,
    float *yheat, float *ypme, float *ytaux, float *ytauy, float *yust3,
    float *yicesnow, float *ycliced);

void ocean_mksst_begin(
    int32_t nhor, float dtmix, float crhos, float cps,
    float tfreeze, const float *yls, const float *yiced,
    const float *yheat, const float *yfldo, const float *ymld,
    const float *ysst, double *zsst, float *yiflux);

void ocean_mksst_finish(
    int32_t nhor, int32_t nocean, int32_t nfluko,
    float dtmix, float crhos, float cps, float tfreeze,
    const float *yls, const float *yiced, const float *ymld,
    double *zsst, float *yiflux, float *ysst);

int32_t ocean_finish_step(
    int32_t nhor, int32_t nocean, int32_t nfluko, int32_t nstep,
    int32_t nout, float dtmix, float taunc, float crhos, float crhoi,
    float cps, float clfi, float tfreeze, const float *yls,
    const float *yclsst2, const float *ymld, const float *yiced,
    const float *ycliced, const float *yheat, const float *yfldo,
    float *ysst, float *yiflux, float *yifluxr,
    float *yheata, float *yifluxa, float *yfldoa, int32_t *naccuout);

void ocean_reset_diagnostics(
    int32_t nhor, float *yheata, float *yifluxa,
    float *yfldoa, int32_t *naccuout);

#endif
