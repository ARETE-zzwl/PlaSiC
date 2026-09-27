#ifndef PLASIC_GAUSSMOD_KERNELS_H
#define PLASIC_GAUSSMOD_KERNELS_H

#include <stdint.h>

double legendre_norm(int32_t degree, double point);
void gauss_inigau(
    int32_t latitudes, double *abscissas, double *weights);

#endif
