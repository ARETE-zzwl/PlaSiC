#include "p_earth_kernels.h"

/*
 * 设置 PlaSiC 的地球默认参数。每个输出指针只代表一种常数，名字同时给出
 * 物理意义或单位：周期以 seconds 表示，行星半径以 m 表示。三个
 * saturation_vapor_pressure_* 参数共同进入模式的饱和水汽压经验式；
 * kappa=Rd/cp，无量纲。函数只写常数，不借用任何参数保存中间结果。
 */
void planet_defaults(
    float *sidereal_day_seconds, float *solar_day_seconds, float *kappa,
    float *dry_air_gas_constant,
    float *saturation_vapor_pressure_reference,
    float *saturation_vapor_pressure_exponent,
    float *saturation_vapor_pressure_temperature_offset,
    float *robert_asselin_filter_weight, float *gravity,
    float *planet_radius)
{
    /*
     * 地球恒星日长度，单位为秒。
     * 恒星日略短于 24 小时。
     */
    *sidereal_day_seconds = 86164.0916f;

    /* 地球平均太阳日长度，单位为秒，即 24 小时。 */
    *solar_day_seconds = 86400.0f;

    /*
     * 设置大气热力学及水汽相关参数。
     */

    /*
     * 干空气热力学常数 κ，近似等于 R/cp。
     * 该参数常用于位温、绝热过程和垂直坐标变换。
     */
    *kappa = 0.286f;

    /* 干空气比气体常数，单位为 J/(kg·K)。 */
    *dry_air_gas_constant = 287.0f;

    /*
     * 饱和水汽压公式中的基准系数。
     * 610.78 Pa 接近水在 0°C 附近的饱和水汽压。
     */
    *saturation_vapor_pressure_reference = 610.78f;

    *saturation_vapor_pressure_exponent = 17.2693882f;

    *saturation_vapor_pressure_temperature_offset = 35.86f;

    *robert_asselin_filter_weight = 0.1f;

    /* 标准重力加速度，单位为 m/s²。 */
    *gravity = 9.80665f;

    /*
     * PlaSiC 使用的地球等效半径，单位为 m。
     * 该值可能作为球形地球近似下的模型半径。
    */
    *planet_radius = 6371220.0f;
}
