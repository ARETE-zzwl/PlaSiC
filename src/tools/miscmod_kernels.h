#ifndef PLASIC_MISCMOD_KERNELS_H
#define PLASIC_MISCMOD_KERNELS_H

#include <stdint.h>

void misc_initialize(
    int32_t longitude_count, int32_t latitude_count,
    const double *gaussian_weights, double two_pi, float rotation_rate,
    float *nudging_time, float *grid_weights);

void misc_prepare_fixer(
    int32_t horizontal_count, int32_t level_count,
    const float *humidity, const float *sigma_thickness,
    const float *surface_pressure, const float *grid_weights,
    float *corrected_humidity, float *moisture_needed,
    float *moisture_available);

float misc_sum(int32_t count, const float *values);

void misc_fix_columns(
    int32_t horizontal_count, int32_t level_count,
    float *corrected_humidity, float *moisture_needed,
    float *moisture_available);

void misc_fix_latitudes(
    int32_t longitude_count, int32_t latitude_count, int32_t level_count,
    float *corrected_humidity, float *moisture_needed,
    float *moisture_available);

void misc_local_totals(
    int32_t horizontal_count, const float *moisture_needed,
    const float *moisture_available, float totals[2]);

void misc_finish_fixer(
    int32_t horizontal_count, int32_t level_count,
    float global_moisture_needed, float global_moisture_available,
    const float *surface_pressure, float twice_timestep,
    float reference_pressure, float *corrected_humidity,
    float *humidity, float *humidity_tendency);

#endif
