#ifndef PLASIC_P_EARTH_KERNELS_H
#define PLASIC_P_EARTH_KERNELS_H

/* 在读取配置前写入地球自转和热力学默认参数。 */
void planet_defaults(
    float *sidereal_day_seconds, float *solar_day_seconds, float *kappa,
    float *dry_air_gas_constant,
    float *saturation_vapor_pressure_reference,
    float *saturation_vapor_pressure_exponent,
    float *saturation_vapor_pressure_temperature_offset,
    float *robert_asselin_filter_weight, float *gravity,
    float *planet_radius);

#endif
