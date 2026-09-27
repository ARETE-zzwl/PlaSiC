#ifndef PLASIC_SEAMOD_KERNELS_H
#define PLASIC_SEAMOD_KERNELS_H

#include <stdint.h>

/*
 * C numerical kernels for seamod.f90.
 *
 * All arrays are contiguous float storage. Two-dimensional fields retain the
 * model's column-major layout: (horizontal, level) is stored at horizontal +
 * level * nhor. Level indices passed here are zero based.
 */

void sea_initialize(
    int32_t nhor, int32_t surface_level, int32_t nrestart,
    float rdbrv, float ra1, float ra2, float ra4, float tmelt, float psurf,
    float albsea, float albice, float dz0sea, float dz0ice,
    float dwetfacsea, float dwetfacice,
    const float *dls, const float *cts, const float *csst,
    const float *cmld, const float *cicec, const float *ciced,
    const float *csnow, float *dts, float *dicec, float *diced,
    float *dsnow, float *dsst, float *dmld, float *dt, float *dq,
    float *dqs, float *dwetfac, float *dalb, float *dz0);

int32_t sea_prepare_step(
    int32_t nhor, int32_t surface_level, int32_t nstep,
    int32_t ncpl_atmos_ice, int32_t nkits, int32_t *naccua,
    const float *dls, const float *dshfl, const float *dswfl,
    const float *dlwfl, const float *dlhfl, const float *dprl,
    const float *dprc, const float *devap, const float *dprs,
    const float *dtaux, const float *dtauy, const float *dust3,
    const float *dshdt, const float *dlhdt, float *cheata, float *cpmea,
    float *cprsa, float *ctauxa, float *ctauya, float *cust3a, float *cshfla,
    float *cshdta, float *clhfla, float *clhdta, float *cswfla,
    float *clwfla);

void sea_finish_coupling(
    int32_t nhor, int32_t *naccua, const float *dls,
    const float *cts, const float *csst, const float *cmld,
    const float *cicec, const float *ciced, const float *csnow,
    const float *csmelt, const float *csndch,
    float *dts, float *dsst, float *dmld, float *dicec, float *diced,
    float *dsnow, float *dsmelt, float *dsndch,
    float *cheata, float *cpmea, float *cprsa, float *ctauxa,
    float *ctauya, float *cust3a, float *cshfla,
    float *cshdta, float *clhfla, float *clhdta, float *cswfla,
    float *clwfla);

void sea_update_surface(
    int32_t nhor, int32_t surface_level,
    float rdbrv, float ra1, float ra2, float ra4, float tmelt,
    float charnock, float gascon, float gravity,
    float albsea, float albice, float dz0sea, float dz0ice,
    float dwetfacsea, float dwetfacice,
    const float *dls, const float *dts, const float *dicec,
    float *dt, float *dq, const float *dp,
    const float *dtaux, const float *dtauy,
    float *dqs, float *dwetfac, float *dalb, float *dz0);

#endif
