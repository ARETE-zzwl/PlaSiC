#ifndef PLASIC_SURFMOD_KERNELS_H
#define PLASIC_SURFMOD_KERNELS_H

#include <stdint.h>

/* Expands one or twelve monthly fields to the model's 14-month layout. */
void surface_expand_months(
    int32_t point_count, int32_t months_read, float *field);

/* Non-dimensionalises spectral orography in place (divide by cv^2). */
void surface_prepare_spectral(
    int32_t coefficient_count, float cv, float *spectral_orography);

#endif
