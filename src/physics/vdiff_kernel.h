#ifndef PLASIC_VDIFF_KERNEL_H
#define PLASIC_VDIFF_KERNEL_H

/*
 * C implementation of the atmospheric vertical-diffusion kernel.
 *
 * All two-dimensional fields retain the model's column-major layout:
 * field(horizontal_index, level_index).  Therefore adjacent horizontal
 * points are contiguous and the C offset is horizontal + level * nhor.
 */
void vdiff(
    int nhor,
    int nlev,
    int ndheat,
    float zumin,
    float vdiff_lamm,
    float vdiff_b,
    float vdiff_c,
    float vdiff_d,
    float deltsec2,
    float ga,
    float gascon,
    float kap,
    float rdbrv,
    float cpd,
    float cpv_cpd_minus1,
    const float *sigma,
    const float *sigmah,
    const float *dsigma,
    const float *dt,
    const float *dq,
    const float *du,
    const float *dv,
    float *dtdt,
    float *dqdt,
    float *dudt,
    float *dvdt);

#endif
