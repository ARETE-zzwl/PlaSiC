#ifndef PLASIC_RAINMOD_KERNELS_H
#define PLASIC_RAINMOD_KERNELS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void rain_initialize_defaults(
    int32_t truncation, int32_t nlev, const float *sigma,
    int32_t *shallow_convection, float *critical_humidity);

float rain_cloud_factor(float first_threshold, float second_threshold);

void rain_zero_precipitation(
    int32_t nhor, float *dprl, float *dprc, float *dprs);

void rain_scale_dqt(
    int32_t nhor, int32_t nlev, float psurf, float omega,
    const float *surface_pressure, float *dqt);

void rain_lsp_prepare(
    int32_t nhor, float sigma_level, float deltsec2, float lv,
    float cpd, float cpv_cpd_minus1, float rdbrv, float ra1, float ra2,
    float tmelt, float ra4, const float *dq, const float *dqdt,
    const float *dt, const float *dtdt, const float *surface_pressure,
    float *predicted_q, float *predicted_t, float *latent_over_cp,
    float *qsat, float *correction, float *effective_qsat);

void rain_lsp_finish(
    int32_t nhor, float sigma_thickness, float sigma_level,
    float deltsec2, float rdbrv, float ra1, float ra2, float tmelt,
    float ra4, const float *surface_pressure, const float *predicted_q,
    const float *predicted_t, const float *latent_over_cp,
    float *qsat, float *correction, const float *effective_qsat,
    float *dqdt, float *dtdt, float *layer_rain);

int32_t rain_clouds(
    int32_t nhor, int32_t nlev, float solar_day, float deltsec2,
    float rdbrv, float ra1, float ra2, float tmelt, float ra4,
    float clwfac, float clwcrit1, float clwcrit2, float gascon,
    float gravity, const int32_t *convective_layer_count,
    const int32_t *convective_layer_flag, const float *convective_rain,
    const float *surface_pressure, const float *sigma,
    const float *sigma_half, const float *sigma_thickness,
    const float *critical_humidity, const float *dt, const float *dtdt,
    const float *dq, const float *dqdt, const float *vertical_velocity,
    float *cloud_cover, float *cloud_liquid_water,
    float *integrated_humidity);

int32_t rain_dry_adjustment(
    int32_t nhor, int32_t nlev, float kap, float deltsec2,
    const float *sigma, const float *sigma_thickness, const float *dt,
    const float *dq, float *dtdt, float *dqdt);

int32_t rain_shallow_convection(
    int32_t nhor, int32_t nlev, int32_t mix_momentum,
    int32_t dissipative_heating, float shallow_diffusion,
    float shallow_top_pressure, float gravity, float deltsec2,
    float gascon, float kap, float rdbrv, float ra1, float ra2,
    float tmelt, float ra4, float cpd, float cpv_cpd_minus1, const float *sigma,
    const float *sigma_half, const float *sigma_thickness,
    const float *surface_pressure, const float *du, const float *dv,
    const float *dq, float *dudt, float *dvdt, const float *temperature,
    const float *humidity, const int32_t *shallow_flag,
    const int32_t *lift_level, int32_t *top_level,
    float *temperature_tendency, float *humidity_tendency);

int32_t rain_kuo(
    int32_t nhor, int32_t nlev, int32_t beta_mode,
    int32_t surface_convection, int32_t shallow_convection,
    int32_t mix_momentum, int32_t beta_exponent,
    int32_t dissipative_heating, float prescribed_beta,
    float critical_relative_humidity, float deep_pressure,
    float shallow_top_pressure, float shallow_diffusion, float deltsec2,
    float kap, float rdbrv, float ra1, float ra2, float tmelt,
    float ra4, float ls, float lv, float cpd, float cpv_cpd_minus1,
    float gravity, float gascon, const float *sigma,
    const float *sigma_half, const float *sigma_thickness,
    const float *surface_pressure, const float *dt, const float *dq,
    const float *du, const float *dv, const float *dqt, float *dtdt,
    float *dqdt, float *dudt, float *dvdt,
    int32_t *convective_layer_flag, int32_t *convective_layer_count,
    float *layer_convective_rain, float *layer_convective_snow);

int32_t rain_fall(
    int32_t nhor, int32_t nlev, int32_t evaporate_precipitation,
    float gamma, float deltsec2, float ls, float lv, float cpd,
    float cpv_cpd_minus1, float tmelt, float rdbrv, float ra1, float ra2,
    float ra4, float gravity, const float *sigma,
    const float *sigma_thickness, const float *surface_pressure,
    const float *dt, const float *dq, const float *layer_convective_rain,
    const float *layer_large_scale_rain,
    const float *layer_convective_snow,
    const float *layer_large_scale_snow, float *dtdt, float *dqdt,
    float *surface_convective_rain, float *surface_large_scale_rain,
    float *surface_snow);

#ifdef __cplusplus
}
#endif

#endif
