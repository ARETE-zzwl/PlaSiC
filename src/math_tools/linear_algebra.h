#ifndef PLASIC_LINEAR_ALGEBRA_H
#define PLASIC_LINEAR_ALGEBRA_H

#include <stdint.h>

/* Calculate the inverse of a square single-precision matrix. */
void matrix_inverse(
    int32_t size, const float *matrix, float *inverse_matrix);

/* Solve a symmetric positive-definite tridiagonal system in place. */
int32_t solve_spd_tridiagonal(
    int32_t size, float *diagonal, float *off_diagonal,
    float *right_hand_side);

#endif
