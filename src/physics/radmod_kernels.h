#ifndef PLASIC_RADMOD_KERNELS_H
#define PLASIC_RADMOD_KERNELS_H

#include <stdint.h>

void rad_zero_fluxes(
    int32_t nhor, int32_t nlep, float *dfu, float *dfd, float *dftu,
    float *dftd, float *dswfl, float *dlwfl);

void rad_compute_tendencies(
    int32_t nhor, int32_t nlev, float gravity, float dry_air_heat_capacity,
    float cpv_cpd_minus1, const float *dsigma, const float *surface_pressure,
    const float *humidity, const float *shortwave_flux,
    const float *longwave_flux, float *total_flux,
    float *temperature_tendency, float *shortwave_tendency,
    float *longwave_tendency, float *radiative_tendency);

void rad_solar_angles(
    int32_t nlon, int32_t nlat, int32_t minute_of_day,
    float pi, float two_pi, float declination, float distance_factor,
    const double *sin_latitude, const float *cos_latitude,
    float *cos_zenith, float *stored_distance_factor);

int32_t rad_synthetic_ozone(
    int32_t nlon, int32_t nlat, int32_t nlev, float calendar_day,
    int32_t days_per_year, float two_pi,
    float a0, float a1, float seasonal_amplitude, float profile_height,
    float profile_thickness, float seasonal_offset,
    float gas_constant, float gravity, const double *sin_latitude,
    const float *sigma_half, const float *dsigma,
    const float *surface_pressure, const float *temperature,
    float *ozone);

int32_t rad_shortwave(
    int32_t nhor, int32_t nlev, int32_t rayleigh_scattering,
    float gravity, float solar_constant, float distance_factor,
    float cloud_tuning_1, float cloud_tuning_2, float cloud_tuning_3,
    const float *sigma, const float *dsigma,
    const float *surface_pressure, float *surface_albedo,
    const float *land_mask, const float *sea_ice,
    const float *humidity, const float *cloud_liquid,
    const float *cloud_cover, const float *ozone,
    const float *cos_zenith, const float *temperature,
    float *upward_flux, float *downward_flux, float *net_flux);

int32_t rad_longwave(
    int32_t nhor, int32_t nlev, float gravity,
    float water_continuum, float cloud_absorption,
    const float *sigma, const float *sigma_half,
    const float *dsigma, const float *surface_pressure,
    const float *land_mask, const float *humidity,
    const float *temperature, const float *cloud_cover,
    const float *cloud_liquid, const float *ozone,
    const float *carbon_dioxide, float *upward_flux,
    float *downward_flux, float *net_flux);

#endif
