#include "surfmod_kernels.h"

#include <stddef.h>
#include <string.h>

void surface_expand_months(
    int32_t point_count, int32_t months_read, float *field)
{
    int32_t month;

    if (months_read == 12)
    {
        for (month = 11; month >= 0; --month)
        {
            memmove(
                field + (size_t)(month + 1) * (size_t)point_count,
                field + (size_t)month * (size_t)point_count,
                (size_t)point_count * sizeof(*field));
        }
        memcpy(
            field + (size_t)13 * (size_t)point_count,
            field + (size_t)point_count,
            (size_t)point_count * sizeof(*field));
        memcpy(
            field,
            field + (size_t)12 * (size_t)point_count,
            (size_t)point_count * sizeof(*field));
    }
    else if (months_read == 1)
    {
        for (month = 1; month < 14; ++month)
        {
            memcpy(
                field + (size_t)month * (size_t)point_count,
                field,
                (size_t)point_count * sizeof(*field));
        }
    }
}

void surface_prepare_spectral(
    int32_t coefficient_count, float cv, float *spectral_orography)
{
    const float cv_squared = cv * cv;
    int32_t coefficient;

    for (coefficient = 0; coefficient < coefficient_count; ++coefficient)
    {
        spectral_orography[coefficient] /= cv_squared;
    }
}
