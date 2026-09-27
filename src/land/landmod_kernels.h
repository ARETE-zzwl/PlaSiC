#ifndef PLASIC_LANDMOD_KERNELS_H
#define PLASIC_LANDMOD_KERNELS_H

#include <stdint.h>

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
    float *atmosphere_surface_humidity);

void land_apply_glacier_surface(
    int32_t point_count, int32_t initialize_snow, float maximum_snow_depth,
    float glacier_albedo_minimum, float glacier_albedo_maximum,
    float melting_temperature, const float *land_mask,
    const float *glacier_mask, const float *surface_temperature,
    float *snow_depth, float *diagnostic_snow, float *surface_albedo,
    float *surface_wetness);

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
    float *soil_heat_flux);

int32_t land_mktsoil(
    int32_t point_count, int32_t soil_levels, float timestep,
    float melting_temperature, const float *land_mask,
    const float *glacier_mask, const float *surface_heat_flux,
    const float *soil_layer_depth, const float *heat_capacity,
    const float *conductivity, float *soil_temperature);

void land_update_soil_water(
    int32_t point_count, float timestep, const float *land_mask,
    const float *soil_water_flux, const float *field_capacity,
    float *soil_water);

void land_prescribed_temperature(
    int32_t point_count, const float *land_mask, const float *evaporation,
    const float *large_scale_precipitation,
    const float *convective_precipitation,
    const float *temperature_climatology, float *soil_water_flux,
    float *surface_temperature);

void land_prescribed_water(
    int32_t point_count, const float *land_mask,
    const float *wetness_climatology, float *soil_water);

void land_interpolate_climatology(
    int32_t point_count, int32_t first_month, int32_t second_month,
    float second_weight, const float *temperature_climatology,
    const float *wetness_climatology, float *interpolated_temperature,
    float *interpolated_wetness);

void land_interpolate_albedo(
    int32_t point_count, int32_t first_month, int32_t second_month,
    float second_weight, const float *albedo_climatology,
    float *interpolated_albedo);

#endif
