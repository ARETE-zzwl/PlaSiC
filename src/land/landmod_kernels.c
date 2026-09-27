#include "landmod_kernels.h"

#include "math_tools/linear_algebra.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

static size_t land_index(int32_t point, int32_t level, int32_t point_count)
{
    return (size_t)level * (size_t)point_count + (size_t)point;
}

static float land_min(float left, float right)
{
    return left < right ? left : right;
}

static float land_max(float left, float right)
{
    return left > right ? left : right;
}

void land_surface_exchange(
    int32_t point_count, int32_t soil_levels, int32_t initialize_roughness,
    int32_t use_constant_pressure, float constant_pressure, float rd_over_rv,
    float saturation_a1, float saturation_a2, float saturation_a4,
    float melting_temperature, float snow_albedo_minimum,
    float snow_albedo_maximum, float full_wetness_fraction,
    const float *pressure, const float *land_mask,
    const float *surface_temperature, const float *snow_depth,
    const float *soil_water, const float *field_capacity,
    const float *background_albedo, const float *soil_temperature,
    const float *roughness_climatology, float *previous_surface_temperature,
    float *surface_humidity, float *diagnostic_snow, float *surface_albedo,
    float *surface_wetness, float *roughness, float *soil_level_1,
    float *soil_level_2, float *soil_level_3, float *soil_level_4,
    float *soil_level_5, float *atmosphere_surface_temperature,
    float *atmosphere_surface_humidity)
{
    int32_t point;

    for (point = 0; point < point_count; ++point)
    {
        float humidity;
        float local_pressure;

        if (land_mask[point] <= 0.0f)
        {
            continue;
        }

        previous_surface_temperature[point] = surface_temperature[point];
        local_pressure =
            use_constant_pressure != 0 ? constant_pressure : pressure[point];
        humidity =
            rd_over_rv * saturation_a1 *
            expf(saturation_a2 *
                 (surface_temperature[point] - melting_temperature) /
                 (surface_temperature[point] - saturation_a4)) /
            local_pressure;
        humidity =
            humidity /
            (1.0f - (1.0f / rd_over_rv - 1.0f) * humidity);
        surface_humidity[point] = humidity;
        diagnostic_snow[point] = snow_depth[point];

        if (diagnostic_snow[point] > 0.0f)
        {
            float delta =
                (snow_albedo_maximum - snow_albedo_minimum) *
                (surface_temperature[point] - 263.16f) /
                (melting_temperature - 263.16f);
            float snow_albedo =
                land_max(
                    snow_albedo_minimum,
                    land_min(
                        snow_albedo_maximum,
                        snow_albedo_maximum - delta));
            surface_albedo[point] =
                background_albedo[point] +
                (snow_albedo - background_albedo[point]) *
                    diagnostic_snow[point] /
                    (diagnostic_snow[point] + 0.01f);
            surface_wetness[point] = 1.0f;
        }
        else
        {
            surface_albedo[point] = background_albedo[point];
            if (field_capacity[point] > 0.0f)
            {
                surface_wetness[point] =
                    land_min(1.0f, soil_water[point] /
                                       (full_wetness_fraction *
                                        field_capacity[point]));
            }
        }

        if (initialize_roughness != 0)
        {
            roughness[point] = roughness_climatology[point];
        }

        soil_level_1[point] =
            soil_temperature[land_index(point, 0, point_count)];
        if (soil_levels > 2)
        {
            soil_level_2[point] =
                soil_temperature[land_index(point, 1, point_count)];
        }
        if (soil_levels > 3)
        {
            soil_level_3[point] =
                soil_temperature[land_index(point, 2, point_count)];
        }
        if (soil_levels > 4)
        {
            soil_level_4[point] =
                soil_temperature[land_index(point, 3, point_count)];
        }
        soil_level_5[point] =
            soil_temperature[land_index(point, soil_levels - 1, point_count)];
        atmosphere_surface_temperature[point] = surface_temperature[point];
        atmosphere_surface_humidity[point] = humidity;
    }
}

void land_apply_glacier_surface(
    int32_t point_count, int32_t initialize_snow, float maximum_snow_depth,
    float glacier_albedo_minimum, float glacier_albedo_maximum,
    float melting_temperature, const float *land_mask,
    const float *glacier_mask, const float *surface_temperature,
    float *snow_depth, float *diagnostic_snow, float *surface_albedo,
    float *surface_wetness)
{
    int32_t point;

    for (point = 0; point < point_count; ++point)
    {
        if (glacier_mask[point] > 0.5f && land_mask[point] > 0.0f)
        {
            float delta;
            if (initialize_snow != 0)
            {
                snow_depth[point] = land_max(maximum_snow_depth, 0.0f);
                diagnostic_snow[point] = snow_depth[point];
            }
            delta =
                (glacier_albedo_maximum - glacier_albedo_minimum) *
                (surface_temperature[point] - 263.16f) /
                (melting_temperature - 263.16f);
            surface_albedo[point] =
                land_max(glacier_albedo_minimum,
                         land_min(glacier_albedo_maximum,
                                  glacier_albedo_maximum - delta));
            surface_wetness[point] = 1.0f;
        }
    }
}

int32_t land_mktsoil(
    int32_t point_count, int32_t soil_levels, float timestep,
    float melting_temperature, const float *land_mask,
    const float *glacier_mask, const float *surface_heat_flux,
    const float *soil_layer_depth, const float *heat_capacity,
    const float *conductivity, float *soil_temperature)
{
    float *memory;
    float *diagonal;
    float *off_diagonal;
    float *right_hand_side;
    int32_t point;
    int32_t level;

    if (point_count <= 0 || soil_levels < 2)
    {
        return 1;
    }
    memory = (float *)malloc((size_t)soil_levels * 3U * sizeof(float));
    if (memory == NULL)
    {
        return 1;
    }
    diagonal = memory;
    off_diagonal = diagonal + soil_levels;
    right_hand_side = off_diagonal + soil_levels;

    for (point = 0; point < point_count; ++point)
    {
        int32_t status;

        if (land_mask[point] <= 0.0f)
        {
            continue;
        }

        /* Conductance between adjacent soil-layer centres. */
        for (level = 0; level < soil_levels - 1; ++level)
        {
            size_t current = land_index(point, level, point_count);
            size_t below = land_index(point, level + 1, point_count);
            float interface_conductivity =
                2.0f * conductivity[current] * conductivity[below] /
                (conductivity[current] * soil_layer_depth[below] +
                 conductivity[below] * soil_layer_depth[current]);

            off_diagonal[level] = -interface_conductivity;
        }

        /* Assemble A T(new) = b for this vertical soil column. */
        for (level = 0; level < soil_levels; ++level)
        {
            size_t index = land_index(point, level, point_count);
            float layer_capacity =
                heat_capacity[index] * soil_layer_depth[index] / timestep;

            diagonal[level] = layer_capacity;
            if (level > 0)
            {
                diagonal[level] -= off_diagonal[level - 1];
            }
            if (level < soil_levels - 1)
            {
                diagonal[level] -= off_diagonal[level];
            }
            right_hand_side[level] =
                layer_capacity * soil_temperature[index];
        }
        right_hand_side[0] += surface_heat_flux[point];

        status = solve_spd_tridiagonal(
            soil_levels, diagonal, off_diagonal, right_hand_side);
        if (status != 0)
        {
            free(memory);
            return status;
        }

        for (level = 0; level < soil_levels; ++level)
        {
            size_t index = land_index(point, level, point_count);
            soil_temperature[index] = right_hand_side[level];
            if (glacier_mask[point] > 0.5f)
            {
                soil_temperature[index] =
                    land_min(soil_temperature[index], melting_temperature);
            }
        }
    }

    free(memory);
    return 0;
}

int32_t land_tands(
    int32_t point_count, int32_t soil_levels, float timestep,
    float top_layer_depth, float maximum_snow_depth,
    float melting_temperature, float latent_heat_sublimation,
    float latent_heat_vaporization, float snow_density,
    float soil_conductivity, float ice_conductivity,
    float snow_conductivity, float soil_heat_capacity,
    float ice_heat_capacity, float snow_heat_capacity,
    const float *soil_layer_depth, const float *land_mask,
    const float *glacier_mask, const float *sensible_heat_flux,
    const float *latent_heat_flux, const float *radiative_surface_flux,
    const float *snow_precipitation, const float *evaporation,
    const float *large_scale_precipitation,
    const float *convective_precipitation,
    const float *previous_surface_temperature, float *soil_temperature,
    float *surface_temperature, float *snow_depth, float *snow_temperature,
    float *soil_water_flux, float *snow_melt, float *snow_depth_change,
    float *atmospheric_heat_flux, float *melt_heat_flux,
    float *soil_heat_flux)
{
    const float maximum_physical_snow_depth = 1.0f;
    size_t field_size;
    size_t vector_size;
    float *memory;
    float *new_snow_depth;
    float *snow_tendency;
    float *capacity;
    float *conductivity;
    float *layer_depth;
    float *top_capacity;
    float *top_conductivity;
    float *top_soil_depth;
    float *mixed_capacity;
    float *top_snow_depth;
    float *mixed_depth;
    float latent_difference;
    int32_t point;
    int32_t level;
    int32_t status;

    if (point_count <= 0 || soil_levels < 2)
    {
        return 1;
    }
    field_size = (size_t)point_count * (size_t)soil_levels;
    vector_size = (size_t)point_count;
    memory = (float *)calloc(field_size * 3U + vector_size * 8U,
                             sizeof(float));
    if (memory == NULL)
    {
        return 1;
    }
    capacity = memory;
    conductivity = capacity + field_size;
    layer_depth = conductivity + field_size;
    new_snow_depth = layer_depth + field_size;
    snow_tendency = new_snow_depth + vector_size;
    top_capacity = snow_tendency + vector_size;
    top_conductivity = top_capacity + vector_size;
    top_soil_depth = top_conductivity + vector_size;
    mixed_capacity = top_soil_depth + vector_size;
    top_snow_depth = mixed_capacity + vector_size;
    mixed_depth = top_snow_depth + vector_size;
    latent_difference = latent_heat_sublimation - latent_heat_vaporization;

    for (point = 0; point < point_count; ++point)
    {
        soil_water_flux[point] = 0.0f;
        snow_melt[point] = 0.0f;
        snow_depth_change[point] = snow_depth[point];
        soil_heat_flux[point] = 0.0f;
        melt_heat_flux[point] = 0.0f;
        mixed_depth[point] = top_layer_depth;
    }

    for (level = 0; level < soil_levels; ++level)
    {
        for (point = 0; point < point_count; ++point)
        {
            size_t index = land_index(point, level, point_count);
            if (land_mask[point] > 0.0f)
            {
                capacity[index] =
                    ice_heat_capacity * glacier_mask[point] +
                    soil_heat_capacity * (1.0f - glacier_mask[point]);
                conductivity[index] =
                    ice_conductivity * glacier_mask[point] +
                    soil_conductivity * (1.0f - glacier_mask[point]);
                layer_depth[index] = soil_layer_depth[level];
            }
        }
    }

    for (point = 0; point < point_count; ++point)
    {
        if (snow_depth[point] == 0.0f)
        {
            snow_temperature[point] = melting_temperature;
        }
        atmospheric_heat_flux[point] =
            sensible_heat_flux[point] + latent_heat_flux[point] +
            radiative_surface_flux[point];

        if (land_mask[point] > 0.0f)
        {
            if (snow_precipitation[point] > 0.0f)
            {
                snow_tendency[point] = snow_precipitation[point];
            }
            if (snow_depth[point] > 0.0f)
            {
                snow_tendency[point] =
                    snow_tendency[point] + evaporation[point];
            }
            new_snow_depth[point] =
                land_max(0.0f, snow_depth[point] +
                                   snow_tendency[point] * timestep);
            snow_tendency[point] =
                (new_snow_depth[point] - snow_depth[point]) / timestep;
            soil_water_flux[point] =
                evaporation[point] + large_scale_precipitation[point] +
                convective_precipitation[point] - snow_tendency[point];
            snow_depth[point] = new_snow_depth[point];
        }
    }

    for (point = 0; point < point_count; ++point)
    {
        size_t top = land_index(point, 0, point_count);
        if (land_mask[point] > 0.0f)
        {
            top_snow_depth[point] =
                land_min(snow_depth[point] * 1000.0f / snow_density,
                         mixed_depth[point]);
            mixed_capacity[point] =
                capacity[top] * snow_heat_capacity * mixed_depth[point] /
                (capacity[top] * top_snow_depth[point] +
                 snow_heat_capacity *
                     (mixed_depth[point] - top_snow_depth[point]));

            new_snow_depth[point] =
                land_min(snow_depth[point] * 1000.0f / snow_density,
                         maximum_physical_snow_depth) -
                top_snow_depth[point];
            top_soil_depth[point] =
                layer_depth[top] + new_snow_depth[point];
            top_conductivity[point] =
                conductivity[top] * snow_conductivity *
                top_soil_depth[point] /
                (snow_conductivity * layer_depth[top] +
                 conductivity[top] * new_snow_depth[point]);
            top_capacity[point] =
                capacity[top] * snow_heat_capacity * top_soil_depth[point] /
                (capacity[top] * new_snow_depth[point] +
                 snow_heat_capacity * layer_depth[top]);

            surface_temperature[point] =
                (mixed_capacity[point] * mixed_depth[point] *
                     previous_surface_temperature[point] / timestep +
                 atmospheric_heat_flux[point] +
                 2.0f * top_conductivity[point] *
                     soil_temperature[top] / top_soil_depth[point]) /
                (mixed_capacity[point] * mixed_depth[point] / timestep +
                 2.0f * top_conductivity[point] /
                     top_soil_depth[point]);
            soil_heat_flux[point] =
                2.0f * top_conductivity[point] / top_soil_depth[point] *
                (surface_temperature[point] - soil_temperature[top]);
        }
    }

    for (point = 0; point < point_count; ++point)
    {
        size_t top = land_index(point, 0, point_count);
        if (land_mask[point] > 0.0f && snow_depth[point] > 0.0f &&
            surface_temperature[point] > melting_temperature)
        {
            surface_temperature[point] = melting_temperature;
            melt_heat_flux[point] =
                land_max(0.0f,
                         atmospheric_heat_flux[point] -
                             mixed_capacity[point] * mixed_depth[point] *
                                 (surface_temperature[point] -
                                  previous_surface_temperature[point]) /
                                 timestep +
                             2.0f * top_conductivity[point] *
                                 (soil_temperature[top] -
                                  surface_temperature[point]) /
                                 top_soil_depth[point]);
            new_snow_depth[point] =
                land_max(0.0f,
                         snow_depth[point] -
                             melt_heat_flux[point] * timestep /
                                 (latent_difference * 1000.0f));
            snow_tendency[point] =
                (new_snow_depth[point] - snow_depth[point]) / timestep;
            melt_heat_flux[point] =
                -snow_tendency[point] * 1000.0f * latent_difference;
            snow_depth[point] = new_snow_depth[point];
            snow_melt[point] = -snow_tendency[point];
            soil_heat_flux[point] =
                2.0f * top_conductivity[point] / top_soil_depth[point] *
                (surface_temperature[point] - soil_temperature[top]);
            soil_water_flux[point] =
                soil_water_flux[point] - snow_tendency[point];

            top_snow_depth[point] =
                land_min(snow_depth[point] * 1000.0f / snow_density,
                         mixed_depth[point]);
            mixed_capacity[point] =
                capacity[top] * snow_heat_capacity * mixed_depth[point] /
                (capacity[top] * top_snow_depth[point] +
                 snow_heat_capacity *
                     (mixed_depth[point] - top_snow_depth[point]));
            new_snow_depth[point] =
                land_min(snow_depth[point] * 1000.0f / snow_density,
                         maximum_physical_snow_depth) -
                top_snow_depth[point];
            top_soil_depth[point] =
                layer_depth[top] + new_snow_depth[point];
            top_conductivity[point] =
                conductivity[top] * snow_conductivity *
                top_soil_depth[point] /
                (snow_conductivity * layer_depth[top] +
                 conductivity[top] * new_snow_depth[point]);
            top_capacity[point] =
                capacity[top] * snow_heat_capacity * top_soil_depth[point] /
                (capacity[top] * new_snow_depth[point] +
                 snow_heat_capacity * layer_depth[top]);
        }
    }

    for (point = 0; point < point_count; ++point)
    {
        size_t top = land_index(point, 0, point_count);
        if (land_mask[point] > 0.0f && snow_depth[point] == 0.0f &&
            melt_heat_flux[point] > 0.0f)
        {
            surface_temperature[point] =
                (mixed_capacity[point] * mixed_depth[point] *
                     previous_surface_temperature[point] / timestep +
                 atmospheric_heat_flux[point] - melt_heat_flux[point] +
                 2.0f * top_conductivity[point] *
                     soil_temperature[top] / top_soil_depth[point]) /
                (mixed_capacity[point] * mixed_depth[point] / timestep +
                 2.0f * top_conductivity[point] /
                     top_soil_depth[point]);
            soil_heat_flux[point] =
                2.0f * top_conductivity[point] / top_soil_depth[point] *
                (surface_temperature[point] - soil_temperature[top]);
        }
    }

    for (point = 0; point < point_count; ++point)
    {
        size_t top = land_index(point, 0, point_count);
        if (land_mask[point] > 0.0f)
        {
            layer_depth[top] = top_soil_depth[point];
            capacity[top] = top_capacity[point];
            conductivity[top] = top_conductivity[point];
        }
    }

    status = land_mktsoil(
        point_count, soil_levels, timestep, melting_temperature, land_mask,
        glacier_mask, soil_heat_flux, layer_depth, capacity, conductivity,
        soil_temperature);
    if (status != 0)
    {
        free(memory);
        return status;
    }

    for (point = 0; point < point_count; ++point)
    {
        if (land_mask[point] > 0.0f && snow_depth[point] > 0.0f)
        {
            snow_temperature[point] = surface_temperature[point];
        }
    }

    if (maximum_snow_depth > 0.0f)
    {
        for (point = 0; point < point_count; ++point)
        {
            if (land_mask[point] > 0.0f &&
                snow_depth[point] > maximum_snow_depth)
            {
                snow_tendency[point] =
                    (maximum_snow_depth - snow_depth[point]) / timestep;
                snow_depth[point] = maximum_snow_depth;
                soil_water_flux[point] =
                    soil_water_flux[point] - snow_tendency[point];
                surface_temperature[point] =
                    surface_temperature[point] +
                    snow_tendency[point] * 1000.0f * latent_difference *
                        timestep / mixed_capacity[point] /
                        mixed_depth[point];
                snow_melt[point] =
                    snow_melt[point] - snow_tendency[point];
            }
        }
    }

    for (point = 0; point < point_count; ++point)
    {
        if (land_mask[point] > 0.0f)
        {
            snow_depth_change[point] =
                (snow_depth[point] - snow_depth_change[point]) / timestep;
        }
    }

    free(memory);
    return 0;
}

void land_update_soil_water(
    int32_t point_count, float timestep, const float *land_mask,
    const float *soil_water_flux, const float *field_capacity,
    float *soil_water)
{
    int32_t point;

    for (point = 0; point < point_count; ++point)
    {
        if (land_mask[point] > 0.0f)
        {
            soil_water[point] =
                soil_water[point] + timestep * soil_water_flux[point];
            /* Water outside the bucket bounds leaves the modeled system. */
            soil_water[point] =
                land_max(land_min(field_capacity[point], soil_water[point]),
                         0.0f);
        }
    }
}

void land_prescribed_temperature(
    int32_t point_count, const float *land_mask, const float *evaporation,
    const float *large_scale_precipitation,
    const float *convective_precipitation,
    const float *temperature_climatology, float *soil_water_flux,
    float *surface_temperature)
{
    int32_t point;
    for (point = 0; point < point_count; ++point)
    {
        if (land_mask[point] > 0.0f)
        {
            soil_water_flux[point] =
                evaporation[point] + large_scale_precipitation[point] +
                convective_precipitation[point];
            surface_temperature[point] = temperature_climatology[point];
        }
    }
}

void land_prescribed_water(
    int32_t point_count, const float *land_mask,
    const float *wetness_climatology, float *soil_water)
{
    int32_t point;
    for (point = 0; point < point_count; ++point)
    {
        if (land_mask[point] > 0.0f)
        {
            soil_water[point] = wetness_climatology[point];
        }
    }
}

void land_interpolate_climatology(
    int32_t point_count, int32_t first_month, int32_t second_month,
    float second_weight, const float *temperature_climatology,
    const float *wetness_climatology, float *interpolated_temperature,
    float *interpolated_wetness)
{
    float first_weight = 1.0f - second_weight;
    int32_t point;
    for (point = 0; point < point_count; ++point)
    {
        interpolated_temperature[point] =
            first_weight *
                temperature_climatology[land_index(point, first_month, point_count)] +
            second_weight *
                temperature_climatology[land_index(point, second_month, point_count)];
        interpolated_wetness[point] =
            first_weight *
                wetness_climatology[land_index(point, first_month, point_count)] +
            second_weight *
                wetness_climatology[land_index(point, second_month, point_count)];
    }
}

void land_interpolate_albedo(
    int32_t point_count, int32_t first_month, int32_t second_month,
    float second_weight, const float *albedo_climatology,
    float *interpolated_albedo)
{
    float first_weight = 1.0f - second_weight;
    int32_t point;
    for (point = 0; point < point_count; ++point)
    {
        interpolated_albedo[point] =
            first_weight *
                albedo_climatology[land_index(point, first_month, point_count)] +
            second_weight *
                albedo_climatology[land_index(point, second_month, point_count)];
    }
}
