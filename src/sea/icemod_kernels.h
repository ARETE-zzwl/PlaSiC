#ifndef PLASIC_ICEMOD_KERNELS_H
#define PLASIC_ICEMOD_KERNELS_H

#include <stdint.h>

void ice_make_ice(
    int32_t nhor, float tfreeze, float ice_density, float latent_heat_ice,
    float timestep, const float *land_mask, const float *ocean_flux,
    const float *sst, const float *snow_conversion_flux,
    const float *freezing_flux, float *conductive_flux,
    float *ice_melt_flux, float *negative_ice_flux, float *ice_thickness);

void ice_compactness(
    int32_t nhor, const float *old_thickness, const float *new_thickness,
    float *compactness);

void ice_snow(
    int32_t nhor, float minimum_ice_cover, float seawater_density,
    float ice_density, float latent_heat_snow, float timestep,
    const float *land_mask, const float *ice_cover, float *snowfall,
    const float *ice_melt_energy, float *ice_thickness, float *snow,
    float *conductive_flux, float *snow_melt, float *snow_melt_flux,
    float *snow_conversion_flux, float *snow_to_ice);

void ice_conductive_flux(
    int32_t nhor, int32_t ocean_flux_correction, int32_t prognostic_ice,
    float minimum_thickness, float maximum_thickness, float tfreeze,
    float snow_density, float ice_conductivity, float snow_conductivity,
    const float *ice_thickness, const float *snow,
    const float *surface_temperature, const float *sst,
    const float *atmospheric_heat, float *diagnosed_conductive_flux,
    float *conductive_flux, float *freezing_flux);

void ice_skin_temperature(
    int32_t nhor, float minimum_thickness, float maximum_thickness,
    float melting_temperature, float ice_density, float ice_heat_capacity,
    float snow_density, float ice_conductivity, float snow_conductivity,
    float timestep, const float *ice_thickness, const float *snow,
    const float *sst, const float *atmospheric_heat,
    float *surface_temperature, float *surface_storage_flux,
    float *melt_energy);

void ice_interpolate_climatology(
    int32_t nhor, int32_t previous_month1, int32_t previous_month2,
    float previous_weight1, float previous_weight2, int32_t month1,
    int32_t month2, float weight1, float weight2,
    const float *climatological_sst, const float *climatological_cover,
    const float *climatological_thickness, float *previous_sst,
    float *current_sst, float *current_cover, float *current_thickness);

void ice_make_thickness(
    int32_t nlon, int32_t nlat, int32_t nmonths,
    const float *compactness, float *thickness);

#endif
