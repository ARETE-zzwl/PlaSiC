#include "miscmod_kernels.h"

#include <stddef.h>

static float misc_max_zero(float value)
{
    return value > 0.0f ? value : 0.0f;
}
void misc_initialize(
    int32_t longitude_count, int32_t latitude_count,
    const double *gaussian_weights, double two_pi, float rotation_rate,
    float *nudging_time, float *grid_weights)
{
    int32_t latitude;

    for (latitude = 0; latitude < latitude_count; ++latitude)
    {
        const float weight = (float)gaussian_weights[latitude];
        int32_t longitude;

        for (longitude = 0; longitude < longitude_count; ++longitude)
        {
            const size_t point =
                (size_t)latitude * (size_t)longitude_count +
                (size_t)longitude;
            grid_weights[point] = weight;
        }
    }

    *nudging_time = (float)(
        two_pi * (double)(*nudging_time) / (double)rotation_rate);
}

void misc_prepare_fixer(
    int32_t horizontal_count, int32_t level_count,
    const float *humidity, const float *sigma_thickness,
    const float *surface_pressure, const float *grid_weights,
    float *corrected_humidity, float *moisture_needed,
    float *moisture_available)
{
    int32_t level;
    int32_t point;

    for (level = 0; level < level_count; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index =
                (size_t)level * (size_t)horizontal_count + (size_t)point;
            corrected_humidity[index] = humidity[index];
        }
    }

    for (point = 0; point < horizontal_count; ++point)
    {
        float negative_sum = 0.0f;
        float positive_sum = 0.0f;

        for (level = 0; level < level_count; ++level)
        {
            const size_t index =
                (size_t)level * (size_t)horizontal_count + (size_t)point;
            const float value = corrected_humidity[index];
            if (value < 0.0f)
            {
                negative_sum =
                    negative_sum + value * sigma_thickness[level];
            }
            if (value > 0.0f)
            {
                positive_sum =
                    positive_sum + value * sigma_thickness[level];
            }
        }
        moisture_needed[point] =
            -negative_sum * surface_pressure[point] * grid_weights[point];
        moisture_available[point] =
            positive_sum * surface_pressure[point] * grid_weights[point];
    }
}
float misc_sum(int32_t count, const float *values)
{
    float sum = 0.0f;
    int32_t index;

    for (index = 0; index < count; ++index)
    {
        sum = sum + values[index];
    }
    return sum;
}

void misc_fix_columns(
    int32_t horizontal_count, int32_t level_count,
    float *corrected_humidity, float *moisture_needed,
    float *moisture_available)
{
    int32_t point;

    for (point = 0; point < horizontal_count; ++point)
    {
        if (moisture_needed[point] > 0.0f &&
            moisture_available[point] >= moisture_needed[point])
        {
            const float factor =
                (moisture_available[point] - moisture_needed[point]) /
                moisture_available[point];
            int32_t level;

            for (level = 0; level < level_count; ++level)
            {
                const size_t index =
                    (size_t)level * (size_t)horizontal_count +
                    (size_t)point;
                corrected_humidity[index] = misc_max_zero(
                    corrected_humidity[index] * factor);
            }
            moisture_needed[point] = 0.0f;
            moisture_available[point] =
                moisture_available[point] * factor;
        }
    }
}

void misc_fix_latitudes(
    int32_t longitude_count, int32_t latitude_count, int32_t level_count,
    float *corrected_humidity, float *moisture_needed,
    float *moisture_available)
{
    int32_t latitude;

    for (latitude = 0; latitude < latitude_count; ++latitude)
    {
        const int32_t first = latitude * longitude_count;
        const int32_t last = first + longitude_count;
        float negative_sum = 0.0f;
        float positive_sum = 0.0f;
        int32_t point;

        for (point = first; point < last; ++point)
        {
            negative_sum = negative_sum + moisture_needed[point];
        }
        if (negative_sum > 0.0f)
        {
            for (point = first; point < last; ++point)
            {
                positive_sum = positive_sum + moisture_available[point];
            }
            if (positive_sum >= negative_sum)
            {
                const float factor =
                    (positive_sum - negative_sum) / positive_sum;
                int32_t level;

                for (level = 0; level < level_count; ++level)
                {
                    for (point = first; point < last; ++point)
                    {
                        const size_t index =
                            (size_t)level * (size_t)(longitude_count *
                                latitude_count) + (size_t)point;
                        corrected_humidity[index] = misc_max_zero(
                            corrected_humidity[index] * factor);
                    }
                }
                for (point = first; point < last; ++point)
                {
                    moisture_needed[point] = 0.0f;
                    moisture_available[point] =
                        moisture_available[point] * factor;
                }
            }
        }
    }
}

void misc_local_totals(
    int32_t horizontal_count, const float *moisture_needed,
    const float *moisture_available, float totals[2])
{
    totals[0] = misc_sum(horizontal_count, moisture_needed);
    totals[1] = misc_sum(horizontal_count, moisture_available);
}

void misc_finish_fixer(
    int32_t horizontal_count, int32_t level_count,
    float global_moisture_needed, float global_moisture_available,
    const float *surface_pressure, float twice_timestep,
    float reference_pressure, float *corrected_humidity,
    float *humidity, float *humidity_tendency)
{
    int32_t level;

    if (global_moisture_needed > 0.0f)
    {
        const float factor =
            (global_moisture_available - global_moisture_needed) /
            global_moisture_available;
        int32_t point;

        for (level = 0; level < level_count; ++level)
        {
            for (point = 0; point < horizontal_count; ++point)
            {
                const size_t index =
                    (size_t)level * (size_t)horizontal_count +
                    (size_t)point;
                corrected_humidity[index] = misc_max_zero(
                    corrected_humidity[index] * factor);
            }
        }
    }

    for (level = 0; level < level_count; ++level)
    {
        int32_t point;
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index =
                (size_t)level * (size_t)horizontal_count + (size_t)point;
            humidity_tendency[index] = humidity_tendency[index] +
                surface_pressure[point] *
                    (corrected_humidity[index] - humidity[index]) /
                    twice_timestep / reference_pressure;
        }
    }

    for (level = 0; level < level_count; ++level)
    {
        int32_t point;
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index =
                (size_t)level * (size_t)horizontal_count + (size_t)point;
            humidity[index] = corrected_humidity[index];
        }
    }
}
