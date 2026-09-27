#include "rainmod_kernels.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define RAIN_AT(array, horizontal, level, nhor) \
    ((array)[(size_t)(horizontal) + (size_t)(nhor) * (size_t)(level)])

static float rain_min(float left, float right)
{
    return left < right ? left : right;
}

static float rain_max(float left, float right)
{
    return left > right ? left : right;
}

static float rain_saturation_mixing_ratio(
    float temperature, float pressure, float sigma_level, float rdbrv,
    float ra1, float ra2, float tmelt, float ra4)
{
    float value = rdbrv * ra1 *
        expf(ra2 * (temperature - tmelt) / (temperature - ra4)) /
        (pressure * sigma_level);
    return rain_min(value, rdbrv);
}

static float rain_mixing_ratio_correction(float mixing_ratio, float rdbrv)
{
    return 1.0f /
        (1.0f - (1.0f / rdbrv - 1.0f) * mixing_ratio);
}

static float rain_saturation_specific_humidity(
    float temperature, float pressure, float sigma_level, float rdbrv,
    float ra1, float ra2, float tmelt, float ra4)
{
    const float mixing_ratio = rain_saturation_mixing_ratio(
        temperature, pressure, sigma_level, rdbrv, ra1, ra2, tmelt, ra4);
    return mixing_ratio /
        (1.0f - (1.0f / rdbrv - 1.0f) * mixing_ratio);
}

static int rain_adjust_saturated_parcel(
    float pressure, float sigma_level, float rdbrv, float ra1, float ra2,
    float tmelt, float ra4, float ls, float lv, float cpd, float cpv_cpd_minus1,
    float *temperature, float *humidity, float *last_correction)
{
    float mixing_ratio = rain_saturation_mixing_ratio(
        *temperature, pressure, sigma_level, rdbrv, ra1, ra2, tmelt, ra4);
    float correction = rain_mixing_ratio_correction(mixing_ratio, rdbrv);
    float saturation = mixing_ratio * correction;
    float latent_over_cp;
    float saturation_derivative;
    float condensed;

    if (!(saturation < *humidity)) {
        *last_correction = correction;
        return 0;
    }
    if (*temperature < tmelt) {
        latent_over_cp = ls / (cpd * (1.0f + cpv_cpd_minus1 * *humidity));
    } else {
        latent_over_cp = lv / (cpd * (1.0f + cpv_cpd_minus1 * *humidity));
    }
    saturation_derivative = ra2 * (tmelt - ra4) * saturation * correction /
        ((*temperature - ra4) * (*temperature - ra4));
    condensed = (*humidity - saturation) /
        (1.0f + latent_over_cp * saturation_derivative);
    *temperature = *temperature + latent_over_cp * condensed;
    *humidity = *humidity - condensed;

    mixing_ratio = rain_saturation_mixing_ratio(
        *temperature, pressure, sigma_level, rdbrv, ra1, ra2, tmelt, ra4);
    correction = rain_mixing_ratio_correction(mixing_ratio, rdbrv);
    saturation = mixing_ratio * correction;
    saturation_derivative = ra2 * (tmelt - ra4) * saturation * correction /
        ((*temperature - ra4) * (*temperature - ra4));
    condensed = (*humidity - saturation) /
        (1.0f + latent_over_cp * saturation_derivative);
    *temperature = *temperature + latent_over_cp * condensed;
    *humidity = *humidity - condensed;
    *last_correction = correction;
    return 1;
}

static float rain_integer_power(float base, int32_t exponent)
{
    float result = 1.0f;
    int32_t count;

    for (count = 0; count < exponent; ++count) {
        result = result * base;
    }
    return result;
}

void rain_initialize_defaults(
    int32_t truncation, int32_t nlev, const float *sigma,
    int32_t *shallow_convection, float *critical_humidity)
{
    int32_t level;

    if (truncation == 21 && nlev == 5) {
        *shallow_convection = 0;
    }
    for (level = 0; level < nlev; ++level) {
        critical_humidity[level] = rain_max(
            0.85f, rain_max(sigma[level], 1.0f - sigma[level]));
    }
}

float rain_cloud_factor(float first_threshold, float second_threshold)
{
    if (first_threshold > second_threshold) {
        return 1.0f / (first_threshold - second_threshold);
    }
    if (first_threshold < second_threshold) {
        return -1.0f;
    }
    return 1.0f;
}

void rain_zero_precipitation(
    int32_t nhor, float *dprl, float *dprc, float *dprs)
{
    int32_t horizontal;

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        dprl[horizontal] = 0.0f;
        dprc[horizontal] = 0.0f;
        dprs[horizontal] = 0.0f;
    }
}

void rain_scale_dqt(
    int32_t nhor, int32_t nlev, float psurf, float omega,
    const float *surface_pressure, float *dqt)
{
    int32_t level;
    int32_t horizontal;

    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            dqt[index] =
                (dqt[index] * psurf) * omega / surface_pressure[horizontal];
        }
    }
}

void rain_lsp_prepare(
    int32_t nhor, float sigma_level, float deltsec2, float lv,
    float cpd, float cpv_cpd_minus1, float rdbrv, float ra1, float ra2,
    float tmelt, float ra4, const float *dq, const float *dqdt,
    const float *dt, const float *dtdt, const float *surface_pressure,
    float *predicted_q, float *predicted_t, float *latent_over_cp,
    float *qsat, float *correction, float *effective_qsat)
{
    int32_t horizontal;

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        float saturation;

        predicted_q[horizontal] =
            dq[horizontal] + dqdt[horizontal] * deltsec2;
        predicted_t[horizontal] =
            dt[horizontal] + dtdt[horizontal] * deltsec2;
        latent_over_cp[horizontal] =
            lv / (cpd * (1.0f + cpv_cpd_minus1 * dq[horizontal]));
        saturation = rain_saturation_mixing_ratio(
            predicted_t[horizontal], surface_pressure[horizontal],
            sigma_level, rdbrv, ra1, ra2, tmelt, ra4);
        correction[horizontal] =
            rain_mixing_ratio_correction(saturation, rdbrv);
        saturation = saturation * correction[horizontal];
        if (saturation < 0.0f) {
            saturation = 0.0f;
        }
        qsat[horizontal] = saturation;
        effective_qsat[horizontal] = saturation;
    }
}

void rain_lsp_finish(
    int32_t nhor, float sigma_thickness, float sigma_level,
    float deltsec2, float rdbrv, float ra1, float ra2, float tmelt,
    float ra4, const float *surface_pressure, const float *predicted_q,
    const float *predicted_t, const float *latent_over_cp,
    float *qsat, float *correction, const float *effective_qsat,
    float *dqdt, float *dtdt, float *layer_rain)
{
    int32_t horizontal;

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        if (predicted_q[horizontal] > effective_qsat[horizontal]) {
            float local_dq;
            float adjusted_q;
            float adjusted_t;
            float denominator;
            float saturation;

            denominator = 1.0f +
                latent_over_cp[horizontal] * ra2 * (tmelt - ra4) *
                qsat[horizontal] * correction[horizontal] /
                ((predicted_t[horizontal] - ra4) *
                 (predicted_t[horizontal] - ra4));
            local_dq =
                (effective_qsat[horizontal] - predicted_q[horizontal]) /
                denominator;
            adjusted_t = predicted_t[horizontal] -
                local_dq * latent_over_cp[horizontal];
            adjusted_q = predicted_q[horizontal] + local_dq;

            saturation = rain_saturation_mixing_ratio(
                adjusted_t, surface_pressure[horizontal], sigma_level,
                rdbrv, ra1, ra2, tmelt, ra4);
            correction[horizontal] =
                rain_mixing_ratio_correction(saturation, rdbrv);
            saturation = saturation * correction[horizontal];
            qsat[horizontal] = saturation;
            denominator = 1.0f +
                latent_over_cp[horizontal] * ra2 * (tmelt - ra4) *
                saturation * correction[horizontal] /
                ((adjusted_t - ra4) * (adjusted_t - ra4));
            local_dq = (saturation - adjusted_q) / denominator;
            adjusted_t = adjusted_t -
                local_dq * latent_over_cp[horizontal];
            adjusted_q = adjusted_q + local_dq;

            local_dq =
                (adjusted_q - predicted_q[horizontal]) / deltsec2;
            dqdt[horizontal] = dqdt[horizontal] + local_dq;
            dtdt[horizontal] = dtdt[horizontal] +
                (adjusted_t - predicted_t[horizontal]) / deltsec2;
            layer_rain[horizontal] =
                -rain_min(local_dq, 0.0f) * sigma_thickness;
        }
    }
}

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
    float *integrated_humidity)
{
    const float cloud_maximum = 1.0f;
    const float convective_a = 0.245f;
    const float convective_b = 0.125f;
    const float convective_maximum = 0.8f;
    const float convective_minimum = 0.05f;
    const float rain_factor = solar_day * 1000.0f;
    float *total_convective_cloud;
    float *layer_convective_cloud;
    float *height;
    int32_t horizontal;
    int32_t level;

    total_convective_cloud =
        (float *)malloc((size_t)nhor * 3U * sizeof(float));
    height = (float *)calloc(
        (size_t)nhor * (size_t)nlev, sizeof(float));
    if (total_convective_cloud == NULL || height == NULL) {
        free(total_convective_cloud);
        free(height);
        return 1;
    }
    layer_convective_cloud = total_convective_cloud + nhor;

    memset(cloud_liquid_water, 0,
           (size_t)nhor * (size_t)(nlev + 1) * sizeof(float));
    memset(cloud_cover, 0,
           (size_t)nhor * (size_t)nlev * sizeof(float));

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        float total = convective_minimum;
        float layer = convective_minimum;

        if (convective_rain[horizontal] > 0.0f) {
            total = convective_a +
                convective_b * logf(convective_rain[horizontal] * rain_factor);
            total = rain_min(
                convective_maximum, rain_max(convective_minimum, total));
            layer = 1.0f - powf(
                1.0f - total,
                1.0f / (float)convective_layer_count[horizontal]);
        }
        total_convective_cloud[horizontal] = total;
        layer_convective_cloud[horizontal] = layer;
    }

    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const float predicted_t =
                dt[index] + dtdt[index] * deltsec2;
            const float predicted_q =
                dq[index] + dqdt[index] * deltsec2;
            float saturation = rain_saturation_mixing_ratio(
                predicted_t, surface_pressure[horizontal], sigma[level],
                rdbrv, ra1, ra2, tmelt, ra4);
            float relative_humidity;
            float vertical_factor = 1.0f;
            float nonconvective;

            saturation = saturation /
                (1.0f - (1.0f / rdbrv - 1.0f) * saturation);
            relative_humidity = predicted_q / saturation;
            if (convective_layer_flag[index] > 0) {
                cloud_cover[index] = layer_convective_cloud[horizontal];
            }
            if (convective_layer_flag[index] == 2 &&
                total_convective_cloud[horizontal] >= 0.3f) {
                cloud_cover[index] = rain_min(
                    cloud_cover[index] +
                        2.0f * (total_convective_cloud[horizontal] - 0.3f),
                    cloud_maximum);
            }
            relative_humidity =
                relative_humidity * (1.0f - cloud_cover[index]);
            if (sigma[level] > 0.8f && clwfac > 0.0f) {
                vertical_factor = clwfac * rain_max(
                    0.0f, clwcrit1 - vertical_velocity[index]);
                if (clwcrit2 > vertical_velocity[index]) {
                    vertical_factor = 1.0f;
                }
            }
            nonconvective = rain_max(
                0.0f,
                (relative_humidity - critical_humidity[level]) /
                    (1.0f - critical_humidity[level]));
            nonconvective =
                vertical_factor * (nonconvective * nonconvective);
            cloud_cover[index] = rain_min(
                cloud_cover[index] + nonconvective, cloud_maximum);
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        integrated_humidity[horizontal] = 0.0f;
        total_convective_cloud[horizontal] = 0.0f;
    }
    for (level = nlev - 1; level >= 1; --level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const float layer_height =
                -dt[index] * gascon / gravity *
                logf(sigma_half[level - 1] / sigma_half[level]);
            integrated_humidity[horizontal] =
                dq[index] * sigma_thickness[level] +
                integrated_humidity[horizontal];
            height[index] =
                total_convective_cloud[horizontal] + layer_height * 0.5f;
            total_convective_cloud[horizontal] =
                total_convective_cloud[horizontal] + layer_height;
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        const size_t index = (size_t)horizontal;
        const float layer_height =
            -dt[index] * gascon / gravity *
            logf(sigma[0] / sigma_half[0]) * 0.5f;
        integrated_humidity[horizontal] =
            dq[index] * sigma_thickness[0] + integrated_humidity[horizontal];
        height[index] = total_convective_cloud[horizontal] + layer_height;
        cloud_cover[(size_t)horizontal +
                    (size_t)nhor * (size_t)nlev] = 1.0f;
        integrated_humidity[horizontal] =
            integrated_humidity[horizontal] *
            surface_pressure[horizontal] / gravity;
        total_convective_cloud[horizontal] =
            700.0f * logf(1.0f + integrated_humidity[horizontal]);
    }
    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const float scale_height =
                total_convective_cloud[horizontal];
            if (scale_height > 0.0f && cloud_cover[index] > 0.0f) {
                cloud_liquid_water[index] =
                    0.00021f * expf(-height[index] / scale_height) *
                    gascon * dt[index] /
                    (sigma[level] * surface_pressure[horizontal]);
                cloud_liquid_water[index] =
                    rain_max(cloud_liquid_water[index], 1.0e-9f);
            }
            cloud_cover[(size_t)horizontal +
                        (size_t)nhor * (size_t)nlev] *=
                1.0f - cloud_cover[index];
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        const size_t surface =
            (size_t)horizontal + (size_t)nhor * (size_t)nlev;
        cloud_cover[surface] = 1.0f - cloud_cover[surface];
    }

    free(total_convective_cloud);
    free(height);
    return 0;
}

int32_t rain_dry_adjustment(
    int32_t nhor, int32_t nlev, float kap, float deltsec2,
    const float *sigma, const float *sigma_thickness, const float *dt,
    const float *dq, float *dtdt, float *dqdt)
{
    const size_t field_count = (size_t)nhor * (size_t)nlev;
    float *storage = (float *)malloc(
        (field_count * 4U + (size_t)nlev + (size_t)nhor * 8U) *
        sizeof(float));
    int32_t *integer_storage = (int32_t *)malloc(
        (size_t)nhor * 3U * sizeof(int32_t));
    float *zskap;
    float *zt;
    float *zth;
    float *zq;
    float *zqn;
    float *ztht;
    float *zqt;
    float *zsum1;
    float *zsum2;
    float *zsumq1;
    float *zsumq2;
    int32_t *itop;
    int32_t *ibase;
    int32_t *icon;
    int32_t level;
    int32_t horizontal;
    int32_t iteration = 0;

    if (storage == NULL || integer_storage == NULL) {
        free(storage);
        free(integer_storage);
        return 1;
    }
    zt = storage;
    zth = zt + field_count;
    zq = zth + field_count;
    zqn = zq + field_count;
    zskap = zqn + field_count;
    ztht = zskap + nlev;
    zqt = ztht + nhor;
    zsum1 = zqt + nhor;
    zsum2 = zsum1 + nhor;
    zsumq1 = zsum2 + nhor;
    zsumq2 = zsumq1 + nhor;
    itop = integer_storage;
    ibase = itop + nhor;
    icon = ibase + nhor;

    for (level = 0; level < nlev; ++level) {
        zskap[level] = powf(sigma[level], kap);
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            zt[index] = dt[index] + dtdt[index] * deltsec2;
            zth[index] = zt[index] / zskap[level];
            zq[index] = dq[index] + dqdt[index] * deltsec2;
            zqn[index] = zq[index];
        }
    }

    for (;;) {
        int32_t active_count = 0;

        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            icon[horizontal] = 0;
            ibase[horizontal] = 0;
            itop[horizontal] = nlev;
            ztht[horizontal] = 0.0f;
            zqt[horizontal] = 0.0f;
        }
        for (level = 0; level < nlev - 1; ++level) {
            const int32_t next = level + 1;
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                const size_t next_index =
                    (size_t)horizontal + (size_t)nhor * (size_t)next;
                if (zth[index] < zth[next_index]) {
                    ibase[horizontal] = next + 1;
                    itop[horizontal] = level;
                    zsum1[horizontal] =
                        zskap[next] * sigma_thickness[next] +
                        zskap[level] * sigma_thickness[level];
                    zsum2[horizontal] =
                        zskap[next] * sigma_thickness[next] * zth[next_index] +
                        zskap[level] * sigma_thickness[level] * zth[index];
                    zsumq1[horizontal] =
                        sigma_thickness[next] + sigma_thickness[level];
                    zsumq2[horizontal] =
                        sigma_thickness[next] * zqn[next_index] +
                        sigma_thickness[level] * zqn[index];
                    ztht[horizontal] =
                        zsum2[horizontal] / zsum1[horizontal];
                    zqt[horizontal] =
                        zsumq2[horizontal] / zsumq1[horizontal];
                    icon[horizontal] = 1;
                }
            }
        }
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            active_count += icon[horizontal];
        }
        if (active_count == 0) {
            break;
        }
        for (level = nlev - 2; level >= 0; --level) {
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                if (itop[horizontal] == level + 1 &&
                    zth[index] < ztht[horizontal]) {
                    itop[horizontal] = level;
                    zsum1[horizontal] +=
                        zskap[level] * sigma_thickness[level];
                    zsum2[horizontal] +=
                        zskap[level] * sigma_thickness[level] * zth[index];
                    zsumq1[horizontal] += sigma_thickness[level];
                    zsumq2[horizontal] +=
                        sigma_thickness[level] * zqn[index];
                    ztht[horizontal] =
                        zsum2[horizontal] / zsum1[horizontal];
                    zqt[horizontal] =
                        zsumq2[horizontal] / zsumq1[horizontal];
                }
            }
        }
        for (level = 0; level < nlev; ++level) {
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                if (level + 1 > itop[horizontal] &&
                    level + 1 <= ibase[horizontal]) {
                    const size_t index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)level;
                    zth[index] = ztht[horizontal];
                    zqn[index] = zqt[horizontal];
                }
            }
        }
        ++iteration;
        if (iteration > 2 * nlev) {
            free(integer_storage);
            free(storage);
            return 2;
        }
    }

    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            dtdt[index] +=
                (zth[index] * zskap[level] - zt[index]) / deltsec2;
            dqdt[index] += (zqn[index] - zq[index]) / deltsec2;
        }
    }

    free(integer_storage);
    free(storage);
    return 0;
}

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
    float *temperature_tendency, float *humidity_tendency)
{
    const int32_t nlem = nlev - 1;
    const size_t field_count = (size_t)nhor * (size_t)nlev;
    const size_t interface_count = (size_t)nhor * (size_t)nlem;
    float *storage = (float *)calloc(
        field_count * 14U + interface_count + (size_t)nlev * 2U,
        sizeof(float));
    int32_t *minimum_top =
        (int32_t *)malloc((size_t)nhor * sizeof(int32_t));
    float *local_temperature;
    float *local_humidity;
    float *diffusion;
    float *new_temperature;
    float *new_humidity;
    float *local_u;
    float *local_v;
    float *new_u;
    float *new_v;
    float *epsilon;
    float *kinetic_energy;
    float *new_kinetic_energy;
    float *local_temperature_change;
    float *local_u_change;
    float *local_v_change;
    float *sigma_power;
    float *sigma_half_power;
    int32_t level;
    int32_t horizontal;
    float constant_one;
    float constant_two;

    if (storage == NULL || minimum_top == NULL) {
        free(storage);
        free(minimum_top);
        return 1;
    }
    local_temperature = storage;
    local_humidity = local_temperature + field_count;
    diffusion = local_humidity + field_count;
    new_temperature = diffusion + field_count;
    new_humidity = new_temperature + field_count;
    local_u = new_humidity + field_count;
    local_v = local_u + field_count;
    new_u = local_v + field_count;
    new_v = new_u + field_count;
    kinetic_energy = new_v + field_count;
    new_kinetic_energy = kinetic_energy + field_count;
    local_temperature_change = new_kinetic_energy + field_count;
    local_u_change = local_temperature_change + field_count;
    local_v_change = local_u_change + field_count;
    epsilon = local_v_change + field_count;
    sigma_power = epsilon + interface_count;
    sigma_half_power = sigma_power + nlev;

    constant_one = gravity * deltsec2 / gascon;
    constant_two = shallow_diffusion * constant_one * gravity / gascon;
    for (level = 0; level < nlev; ++level) {
        sigma_power[level] = powf(sigma[level], kap);
        sigma_half_power[level] = powf(sigma_half[level], kap);
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            local_temperature[index] = temperature[index];
            local_humidity[index] = humidity[index];
            temperature_tendency[index] = 0.0f;
            humidity_tendency[index] = 0.0f;
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        minimum_top[horizontal] = 1;
    }
    for (level = 1; level < nlem; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0 &&
                surface_pressure[horizontal] * sigma[level] <
                    shallow_top_pressure) {
                minimum_top[horizontal] = level + 1;
            }
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        if (minimum_top[horizontal] > top_level[horizontal]) {
            top_level[horizontal] = minimum_top[horizontal];
        }
    }

    for (level = 0; level < nlem; ++level) {
        const int32_t next = level + 1;
        const float reciprocal_delta_sigma =
            1.0f / (sigma[next] - sigma[level]);
        const float factor = constant_two * sigma_half[level] *
            sigma_half[level] * reciprocal_delta_sigma;
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const size_t next_index =
                (size_t)horizontal + (size_t)nhor * (size_t)next;
            if (shallow_flag[horizontal] > 0 &&
                level + 1 >= top_level[horizontal] &&
                level + 1 < lift_level[horizontal]) {
                const float reciprocal_temperature =
                    (sigma_thickness[next] + sigma_thickness[level]) /
                    (local_temperature[index] * sigma_thickness[next] +
                     local_temperature[next_index] * sigma_thickness[level]);
                diffusion[index] = factor * reciprocal_temperature *
                    reciprocal_temperature;
            }
            if (shallow_flag[horizontal] > 0 &&
                level + 2 == top_level[horizontal]) {
                const float reciprocal_temperature =
                    (sigma_thickness[next] + sigma_thickness[level]) /
                    (local_temperature[index] * sigma_thickness[next] +
                     local_temperature[next_index] * sigma_thickness[level]);
                float lower_saturation = rain_saturation_mixing_ratio(
                    local_temperature[index], surface_pressure[horizontal],
                    sigma[level], rdbrv, ra1, ra2, tmelt, ra4);
                float upper_saturation = rain_saturation_mixing_ratio(
                    local_temperature[next_index],
                    surface_pressure[horizontal], sigma[next], rdbrv, ra1,
                    ra2, tmelt, ra4);
                float correction = rain_mixing_ratio_correction(
                    lower_saturation, rdbrv);
                const float lower_relative_humidity =
                    local_humidity[index] /
                    (lower_saturation * correction);
                float upper_relative_humidity;
                float transition;

                correction = rain_mixing_ratio_correction(
                    upper_saturation, rdbrv);
                upper_relative_humidity =
                    local_humidity[next_index] /
                    (upper_saturation * correction);
                transition = fabsf(
                    (lower_relative_humidity - 0.8f) / 0.2f *
                    (upper_relative_humidity - lower_relative_humidity));
                diffusion[index] = factor * reciprocal_temperature *
                    reciprocal_temperature * transition;
            }
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        if (shallow_flag[horizontal] > 0) {
            const size_t index = (size_t)horizontal;
            epsilon[index] = diffusion[index] /
                (sigma_thickness[0] + diffusion[index]);
            new_humidity[index] = sigma_thickness[0] *
                local_humidity[index] /
                (sigma_thickness[0] + diffusion[index]);
        }
    }
    for (level = 1; level < nlem; ++level) {
        const int32_t previous = level - 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                const size_t previous_index =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)previous;
                const float denominator = sigma_thickness[level] +
                    diffusion[index] + diffusion[previous_index] *
                    (1.0f - epsilon[previous_index]);
                epsilon[index] = diffusion[index] / denominator;
                new_humidity[index] =
                    (local_humidity[index] * sigma_thickness[level] +
                     diffusion[previous_index] *
                         new_humidity[previous_index]) /
                    denominator;
            }
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        if (shallow_flag[horizontal] > 0) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)(nlev - 1);
            const size_t previous_index =
                (size_t)horizontal + (size_t)nhor * (size_t)nlem -
                (size_t)nhor;
            new_humidity[index] =
                (local_humidity[index] * sigma_thickness[nlev - 1] +
                 diffusion[previous_index] * new_humidity[previous_index]) /
                (sigma_thickness[nlev - 1] +
                 diffusion[previous_index] *
                     (1.0f - epsilon[previous_index]));
        }
    }
    for (level = nlem - 1; level >= 0; --level) {
        const int32_t next = level + 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                const size_t next_index =
                    (size_t)horizontal + (size_t)nhor * (size_t)next;
                new_humidity[index] = new_humidity[index] +
                    epsilon[index] * new_humidity[next_index];
            }
        }
    }
    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                humidity_tendency[index] =
                    (new_humidity[index] - local_humidity[index]) / deltsec2;
            }
        }
    }

    for (level = 0; level < nlem; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                diffusion[index] =
                    diffusion[index] * sigma_half_power[level];
            }
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        if (shallow_flag[horizontal] > 0) {
            const size_t index = (size_t)horizontal;
            const float denominator = sigma_thickness[0] +
                diffusion[index] / sigma_power[0];
            epsilon[index] = diffusion[index] / denominator;
            new_temperature[index] =
                sigma_thickness[0] * local_temperature[index] / denominator;
        }
    }
    for (level = 1; level < nlem; ++level) {
        const int32_t previous = level - 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                const size_t previous_index =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)previous;
                const float denominator = sigma_thickness[level] +
                    (diffusion[index] + diffusion[previous_index] *
                        (1.0f - epsilon[previous_index] /
                            sigma_power[previous])) /
                        sigma_power[level];
                epsilon[index] = diffusion[index] / denominator;
                new_temperature[index] =
                    (local_temperature[index] * sigma_thickness[level] +
                     diffusion[previous_index] / sigma_power[previous] *
                         new_temperature[previous_index]) /
                    denominator;
            }
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        if (shallow_flag[horizontal] > 0) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)(nlev - 1);
            const size_t previous_index =
                (size_t)horizontal + (size_t)nhor * (size_t)(nlem - 1);
            new_temperature[index] =
                (local_temperature[index] * sigma_thickness[nlev - 1] +
                 diffusion[previous_index] *
                     new_temperature[previous_index] /
                     sigma_power[nlem - 1]) /
                (sigma_thickness[nlev - 1] +
                 diffusion[previous_index] / sigma_power[nlev - 1] *
                     (1.0f - epsilon[previous_index] /
                         sigma_power[nlem - 1]));
        }
    }
    for (level = nlem - 1; level >= 0; --level) {
        const int32_t next = level + 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                const size_t next_index =
                    (size_t)horizontal + (size_t)nhor * (size_t)next;
                new_temperature[index] = new_temperature[index] +
                    epsilon[index] * new_temperature[next_index] /
                        sigma_power[next];
            }
        }
    }
    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                temperature_tendency[index] =
                    (new_temperature[index] - local_temperature[index]) /
                    deltsec2;
            }
        }
    }

    if (mix_momentum == 1) {
        memcpy(local_u, du, field_count * sizeof(float));
        memcpy(local_v, dv, field_count * sizeof(float));
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index = (size_t)horizontal;
                const float denominator =
                    sigma_thickness[0] + diffusion[index];
                epsilon[index] = diffusion[index] / denominator;
                new_u[index] =
                    sigma_thickness[0] * local_u[index] / denominator;
                new_v[index] =
                    sigma_thickness[0] * local_v[index] / denominator;
            }
        }
        for (level = 1; level < nlem; ++level) {
            const int32_t previous = level - 1;
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                if (shallow_flag[horizontal] > 0) {
                    const size_t index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)level;
                    const size_t previous_index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)previous;
                    const float denominator = sigma_thickness[level] +
                        diffusion[index] + diffusion[previous_index] *
                        (1.0f - epsilon[previous_index]);
                    epsilon[index] = diffusion[index] / denominator;
                    new_u[index] =
                        (local_u[index] * sigma_thickness[level] +
                         diffusion[previous_index] * new_u[previous_index]) /
                        denominator;
                    new_v[index] =
                        (local_v[index] * sigma_thickness[level] +
                         diffusion[previous_index] * new_v[previous_index]) /
                        denominator;
                }
            }
        }
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (shallow_flag[horizontal] > 0) {
                const size_t index =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)(nlev - 1);
                const size_t previous_index =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)(nlem - 1);
                const float denominator = sigma_thickness[nlev - 1] +
                    diffusion[previous_index] *
                    (1.0f - epsilon[previous_index]);
                new_u[index] =
                    (local_u[index] * sigma_thickness[nlev - 1] +
                     diffusion[previous_index] * new_u[previous_index]) /
                    denominator;
                new_v[index] =
                    (local_v[index] * sigma_thickness[nlev - 1] +
                     diffusion[previous_index] * new_v[previous_index]) /
                    denominator;
            }
        }
        for (level = nlem - 1; level >= 0; --level) {
            const int32_t next = level + 1;
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                if (shallow_flag[horizontal] > 0) {
                    const size_t index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)level;
                    const size_t next_index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)next;
                    new_u[index] =
                        new_u[index] + epsilon[index] * new_u[next_index];
                    new_v[index] =
                        new_v[index] + epsilon[index] * new_v[next_index];
                }
            }
        }
        for (level = 0; level < nlev; ++level) {
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                if (shallow_flag[horizontal] > 0) {
                    const size_t index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)level;
                    local_u_change[index] =
                        (new_u[index] - local_u[index]) / deltsec2;
                    local_v_change[index] =
                        (new_v[index] - local_v[index]) / deltsec2;
                    dudt[index] += local_u_change[index];
                    dvdt[index] += local_v_change[index];
                }
            }
        }
        if (dissipative_heating > 0) {
            for (level = 0; level < nlev; ++level) {
                for (horizontal = 0; horizontal < nhor; ++horizontal) {
                    if (shallow_flag[horizontal] > 0) {
                        const size_t index =
                            (size_t)horizontal +
                            (size_t)nhor * (size_t)level;
                        kinetic_energy[index] = 0.5f *
                            (local_u[index] * local_u[index] +
                             local_v[index] * local_v[index]);
                    }
                }
            }
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                if (shallow_flag[horizontal] > 0) {
                    const size_t index = (size_t)horizontal;
                    const float denominator =
                        sigma_thickness[0] + diffusion[index];
                    epsilon[index] = diffusion[index] / denominator;
                    new_kinetic_energy[index] = sigma_thickness[0] *
                        kinetic_energy[index] / denominator;
                }
            }
            for (level = 1; level < nlem; ++level) {
                const int32_t previous = level - 1;
                for (horizontal = 0; horizontal < nhor; ++horizontal) {
                    if (shallow_flag[horizontal] > 0) {
                        const size_t index =
                            (size_t)horizontal +
                            (size_t)nhor * (size_t)level;
                        const size_t previous_index =
                            (size_t)horizontal +
                            (size_t)nhor * (size_t)previous;
                        const float denominator = sigma_thickness[level] +
                            diffusion[index] + diffusion[previous_index] *
                            (1.0f - epsilon[previous_index]);
                        epsilon[index] = diffusion[index] / denominator;
                        new_kinetic_energy[index] =
                            (kinetic_energy[index] *
                                 sigma_thickness[level] +
                             diffusion[previous_index] *
                                 new_kinetic_energy[previous_index]) /
                            denominator;
                    }
                }
            }
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                if (shallow_flag[horizontal] > 0) {
                    const size_t index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)(nlev - 1);
                    const size_t previous_index =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)(nlem - 1);
                    new_kinetic_energy[index] =
                        (kinetic_energy[index] *
                             sigma_thickness[nlev - 1] +
                         diffusion[previous_index] *
                             new_kinetic_energy[previous_index]) /
                        (sigma_thickness[nlev - 1] +
                         diffusion[previous_index] *
                             (1.0f - epsilon[previous_index]));
                }
            }
            for (level = nlem - 1; level >= 0; --level) {
                const int32_t next = level + 1;
                for (horizontal = 0; horizontal < nhor; ++horizontal) {
                    if (shallow_flag[horizontal] > 0) {
                        const size_t index =
                            (size_t)horizontal +
                            (size_t)nhor * (size_t)level;
                        const size_t next_index =
                            (size_t)horizontal +
                            (size_t)nhor * (size_t)next;
                        new_kinetic_energy[index] =
                            new_kinetic_energy[index] +
                            epsilon[index] * new_kinetic_energy[next_index];
                    }
                }
            }
            for (level = 0; level < nlev; ++level) {
                for (horizontal = 0; horizontal < nhor; ++horizontal) {
                    if (shallow_flag[horizontal] > 0) {
                        const size_t index =
                            (size_t)horizontal +
                            (size_t)nhor * (size_t)level;
                        local_temperature_change[index] =
                            -((new_u[index] * new_u[index] -
                               local_u[index] * local_u[index] +
                               new_v[index] * new_v[index] -
                               local_v[index] * local_v[index]) /
                                  deltsec2 -
                              (new_kinetic_energy[index] -
                               kinetic_energy[index]) /
                                  deltsec2) *
                            0.5f / cpd /
                            (1.0f + cpv_cpd_minus1 * dq[index]);
                        temperature_tendency[index] +=
                            local_temperature_change[index];
                        local_temperature_change[index] = 0.0f;
                    }
                }
            }
        }
    }

    free(minimum_top);
    free(storage);
    return 0;
}

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
    float *layer_convective_rain, float *layer_convective_snow)
{
    const int32_t nlem = nlev - 1;
    const size_t field_count = (size_t)nhor * (size_t)nlev;
    float *field_storage = (float *)calloc(field_count * 13U, sizeof(float));
    float *column_storage =
        (float *)calloc((size_t)nhor * 15U, sizeof(float));
    int32_t *integer_storage =
        (int32_t *)malloc((size_t)nhor * 4U * sizeof(int32_t));
    float *environment_q;
    float *environment_t;
    float *moisture_supply;
    float *local_dq;
    float *local_dt;
    float *latent_over_cp;
    float *effective_cp;
    float *shallow_dq;
    float *shallow_dt;
    float *environment_u;
    float *environment_v;
    float *local_du;
    float *local_dv;
    float *parcel_q;
    float *parcel_t;
    float *available_energy;
    float *temperature_energy;
    float *humidity_energy;
    float *beta_sum;
    float *beta_weight;
    float *temperature_factor;
    float *humidity_factor;
    float *top_pressure;
    float *mean_relative_humidity;
    float *momentum_factor;
    float *parcel_u;
    float *parcel_v;
    float *mixed_sigma;
    int32_t *lift_level;
    int32_t *top_level;
    int32_t *shallow_flag;
    int32_t *deep_flag;
    int32_t level;
    int32_t horizontal;
    int32_t status;

    if (field_storage == NULL || column_storage == NULL ||
        integer_storage == NULL) {
        free(field_storage);
        free(column_storage);
        free(integer_storage);
        return 1;
    }
    environment_q = field_storage;
    environment_t = environment_q + field_count;
    moisture_supply = environment_t + field_count;
    local_dq = moisture_supply + field_count;
    local_dt = local_dq + field_count;
    latent_over_cp = local_dt + field_count;
    effective_cp = latent_over_cp + field_count;
    shallow_dq = effective_cp + field_count;
    shallow_dt = shallow_dq + field_count;
    environment_u = shallow_dt + field_count;
    environment_v = environment_u + field_count;
    local_du = environment_v + field_count;
    local_dv = local_du + field_count;

    parcel_q = column_storage;
    parcel_t = parcel_q + nhor;
    available_energy = parcel_t + nhor;
    temperature_energy = available_energy + nhor;
    humidity_energy = temperature_energy + nhor;
    beta_sum = humidity_energy + nhor;
    beta_weight = beta_sum + nhor;
    temperature_factor = beta_weight + nhor;
    humidity_factor = temperature_factor + nhor;
    top_pressure = humidity_factor + nhor;
    mean_relative_humidity = top_pressure + nhor;
    momentum_factor = mean_relative_humidity + nhor;
    parcel_u = momentum_factor + nhor;
    parcel_v = parcel_u + nhor;
    mixed_sigma = parcel_v + nhor;
    lift_level = integer_storage;
    top_level = lift_level + nhor;
    shallow_flag = top_level + nhor;
    deep_flag = shallow_flag + nhor;

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        shallow_flag[horizontal] = 0;
        lift_level[horizontal] = -1;
        top_level[horizontal] = -1;
        convective_layer_count[horizontal] = 0;
        top_pressure[horizontal] = surface_pressure[horizontal];
    }
    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            moisture_supply[index] =
                (dqt[index] + dqdt[index]) * deltsec2;
            environment_q[index] = dq[index] + dqdt[index] * deltsec2;
            environment_t[index] = dt[index] + dtdt[index] * deltsec2;
            latent_over_cp[index] =
                lv / (cpd * (1.0f + cpv_cpd_minus1 * dq[index]));
            effective_cp[index] = cpd * (1.0f + cpv_cpd_minus1 * dq[index]);
            convective_layer_flag[index] = 0;
            layer_convective_rain[index] = 0.0f;
            layer_convective_snow[index] = 0.0f;
        }
    }

    if (surface_convection == 1) {
        const float factor = powf(sigma[nlev - 1], kap);
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            float parcel_temperature =
                RAIN_AT(dt, horizontal, nlev, nhor) * factor;
            const size_t lowest =
                (size_t)horizontal +
                (size_t)nhor * (size_t)(nlev - 1);
            if (moisture_supply[lowest] > 0.0f &&
                parcel_temperature > environment_t[lowest]) {
                const float environment_saturation =
                    rain_saturation_specific_humidity(
                        environment_t[lowest], surface_pressure[horizontal],
                        sigma[nlev - 1], rdbrv, ra1, ra2, tmelt, ra4);
                float parcel_humidity = RAIN_AT(dq, horizontal, nlev, nhor) *
                    rain_min(1.0f,
                             environment_q[lowest] / environment_saturation);
                float correction;

                if (rain_adjust_saturated_parcel(
                        surface_pressure[horizontal], sigma[nlev - 1],
                        rdbrv, ra1, ra2, tmelt, ra4, ls, lv, cpd, cpv_cpd_minus1,
                        &parcel_temperature, &parcel_humidity, &correction)) {
                    const float saturation =
                        rain_saturation_specific_humidity(
                            environment_t[lowest],
                            surface_pressure[horizontal], sigma[nlev - 1],
                            rdbrv, ra1, ra2, tmelt, ra4);
                    lift_level[horizontal] = nlev + 1;
                    top_level[horizontal] = nlev;
                    parcel_t[horizontal] = parcel_temperature;
                    parcel_q[horizontal] = parcel_humidity;
                    local_dq[lowest] =
                        rain_max(saturation - environment_q[lowest], 0.0f);
                    local_dt[lowest] =
                        parcel_temperature - environment_t[lowest];
                    temperature_energy[horizontal] = local_dt[lowest] *
                        sigma_thickness[nlev - 1] * effective_cp[lowest];
                    humidity_energy[horizontal] = local_dq[lowest] *
                        sigma_thickness[nlev - 1] * lv;
                    if (beta_mode > 0) {
                        beta_sum[horizontal] = rain_min(
                            1.0f, environment_q[lowest] / saturation) *
                            sigma_thickness[nlev - 1];
                        beta_weight[horizontal] =
                            sigma_thickness[nlev - 1];
                    }
                    available_energy[horizontal] = moisture_supply[lowest] *
                        sigma_thickness[nlev - 1] * lv;
                }
            }
        }
    }

    for (level = nlem - 1; level >= 0; --level) {
        const int32_t next = level + 1;
        const float factor = powf(sigma[level] / sigma[next], kap);
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const size_t next_index =
                (size_t)horizontal + (size_t)nhor * (size_t)next;
            if (lift_level[horizontal] == -1 &&
                moisture_supply[next_index] > 0.0f) {
                float parcel_temperature = environment_t[next_index] * factor;
                float parcel_humidity = environment_q[next_index];
                float correction;

                if (rain_adjust_saturated_parcel(
                        surface_pressure[horizontal], sigma[level], rdbrv,
                        ra1, ra2, tmelt, ra4, ls, lv, cpd, cpv_cpd_minus1,
                        &parcel_temperature, &parcel_humidity, &correction) &&
                    parcel_temperature > environment_t[index]) {
                    const float saturation =
                        rain_saturation_specific_humidity(
                            environment_t[index], surface_pressure[horizontal],
                            sigma[level], rdbrv, ra1, ra2, tmelt, ra4);
                    float next_temperature;
                    float next_humidity;
                    float next_saturation;

                    lift_level[horizontal] = next + 1;
                    top_level[horizontal] = level + 1;
                    parcel_t[horizontal] = parcel_temperature;
                    parcel_q[horizontal] = parcel_humidity;
                    local_dq[index] =
                        rain_max(saturation - environment_q[index], 0.0f);
                    local_dt[index] =
                        parcel_temperature - environment_t[index];
                    temperature_energy[horizontal] = local_dt[index] *
                        sigma_thickness[level] * effective_cp[index];
                    humidity_energy[horizontal] = local_dq[index] *
                        sigma_thickness[level] * lv;
                    if (beta_mode > 0) {
                        beta_sum[horizontal] = rain_min(
                            1.0f, environment_q[index] / saturation) *
                            sigma_thickness[level];
                        beta_weight[horizontal] = sigma_thickness[level];
                    }

                    next_temperature = environment_t[next_index];
                    next_humidity = environment_q[next_index];
                    next_saturation = rain_saturation_specific_humidity(
                        environment_t[next_index],
                        surface_pressure[horizontal], sigma[next], rdbrv,
                        ra1, ra2, tmelt, ra4);
                    if (environment_q[next_index] > next_saturation) {
                        float latent_over_cp_value;
                        float derivative;
                        float condensed;
                        float mixing_ratio;

                        if (next_temperature < tmelt) {
                            latent_over_cp_value =
                                ls / (cpd * (1.0f + cpv_cpd_minus1 * next_humidity));
                        } else {
                            latent_over_cp_value =
                                lv / (cpd * (1.0f + cpv_cpd_minus1 * next_humidity));
                        }
                        derivative = ra2 * (tmelt - ra4) * next_saturation *
                            correction /
                            ((next_temperature - ra4) *
                             (next_temperature - ra4));
                        condensed = (next_humidity - next_saturation) /
                            (1.0f + latent_over_cp_value * derivative);
                        next_temperature += latent_over_cp_value * condensed;
                        next_humidity -= condensed;
                        mixing_ratio = rain_saturation_mixing_ratio(
                            next_temperature, surface_pressure[horizontal],
                            sigma[next], rdbrv, ra1, ra2, tmelt, ra4);
                        correction = rain_mixing_ratio_correction(
                            mixing_ratio, rdbrv);
                        next_saturation = mixing_ratio * correction;
                        derivative = ra2 * (tmelt - ra4) * next_saturation *
                            correction /
                            ((next_temperature - ra4) *
                             (next_temperature - ra4));
                        condensed = (next_humidity - next_saturation) /
                            (1.0f + latent_over_cp_value * derivative);
                        next_temperature += latent_over_cp_value * condensed;
                        next_humidity -= condensed;
                        local_dq[next_index] = rain_max(
                            next_humidity - environment_q[next_index], 0.0f);
                        local_dt[next_index] =
                            next_temperature - environment_t[next_index];
                        temperature_energy[horizontal] +=
                            local_dt[next_index] * sigma_thickness[next] *
                            effective_cp[next_index];
                        humidity_energy[horizontal] +=
                            local_dq[next_index] * sigma_thickness[next] * lv;
                        if (beta_mode > 0) {
                            beta_sum[horizontal] += rain_min(
                                1.0f,
                                environment_q[next_index] / next_saturation) *
                                sigma_thickness[next];
                            beta_weight[horizontal] += sigma_thickness[next];
                        }
                    }
                    available_energy[horizontal] =
                        moisture_supply[index] * sigma_thickness[level] * lv +
                        moisture_supply[next_index] *
                            sigma_thickness[next] * lv;
                }
            }
        }
    }

    for (level = nlem - 1; level >= 0; --level) {
        const int32_t next = level + 1;
        const float factor = powf(sigma[level] / sigma[next], kap);
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (next + 1 == top_level[horizontal]) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                float parcel_temperature = parcel_t[horizontal] * factor;
                float parcel_humidity = parcel_q[horizontal];
                float correction;

                if (rain_adjust_saturated_parcel(
                        surface_pressure[horizontal], sigma[level], rdbrv,
                        ra1, ra2, tmelt, ra4, ls, lv, cpd, cpv_cpd_minus1,
                        &parcel_temperature, &parcel_humidity, &correction) &&
                    parcel_temperature > environment_t[index]) {
                    const float saturation =
                        rain_saturation_specific_humidity(
                            environment_t[index], surface_pressure[horizontal],
                            sigma[level], rdbrv, ra1, ra2, tmelt, ra4);
                    top_level[horizontal] = level + 1;
                    parcel_t[horizontal] = parcel_temperature;
                    parcel_q[horizontal] = parcel_humidity;
                    local_dq[index] =
                        rain_max(saturation - environment_q[index], 0.0f);
                    local_dt[index] =
                        parcel_temperature - environment_t[index];
                    temperature_energy[horizontal] += local_dt[index] *
                        sigma_thickness[level] * effective_cp[index];
                    humidity_energy[horizontal] += local_dq[index] *
                        sigma_thickness[level] * lv;
                    if (beta_mode > 0) {
                        beta_sum[horizontal] += rain_min(
                            1.0f, environment_q[index] / saturation) *
                            sigma_thickness[level];
                        beta_weight[horizontal] += sigma_thickness[level];
                    }
                    available_energy[horizontal] += moisture_supply[index] *
                        sigma_thickness[level] * lv;
                }
            }
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        available_energy[horizontal] =
            rain_max(0.0f, available_energy[horizontal]);
        if (beta_mode == 0) {
            if (temperature_energy[horizontal] +
                    humidity_energy[horizontal] > 0.0f) {
                temperature_factor[horizontal] =
                    available_energy[horizontal] /
                    (temperature_energy[horizontal] +
                     humidity_energy[horizontal]);
                humidity_factor[horizontal] = temperature_factor[horizontal];
            }
        } else {
            mean_relative_humidity[horizontal] = 0.0f;
            if (beta_weight[horizontal] > 0.0f) {
                mean_relative_humidity[horizontal] =
                    beta_sum[horizontal] / beta_weight[horizontal];
            }
            beta_sum[horizontal] = 1.0f;
            if (mean_relative_humidity[horizontal] >=
                critical_relative_humidity) {
                beta_sum[horizontal] = rain_integer_power(
                    (1.0f - mean_relative_humidity[horizontal]) /
                        (1.0f - critical_relative_humidity),
                    beta_exponent);
            }
            if (beta_mode == 2) {
                beta_sum[horizontal] = prescribed_beta;
            }
            if (temperature_energy[horizontal] > 0.0f) {
                temperature_factor[horizontal] =
                    available_energy[horizontal] *
                    (1.0f - beta_sum[horizontal]) /
                    temperature_energy[horizontal];
            }
            if (humidity_energy[horizontal] > 0.0f) {
                humidity_factor[horizontal] =
                    available_energy[horizontal] * beta_sum[horizontal] /
                    humidity_energy[horizontal];
            }
        }
        if (shallow_convection == 1 && top_level[horizontal] > 0) {
            top_pressure[horizontal] = surface_pressure[horizontal] *
                sigma[top_level[horizontal] - 1];
            if (top_pressure[horizontal] > deep_pressure ||
                available_energy[horizontal] <= 0.0f) {
                shallow_flag[horizontal] = 1;
            }
        }
    }

    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            deep_flag[horizontal] = 0;
            if (lift_level[horizontal] >= level + 1 &&
                top_level[horizontal] <= level + 1 &&
                available_energy[horizontal] > 0.0f &&
                shallow_flag[horizontal] == 0) {
                deep_flag[horizontal] = 1;
            }
            if (deep_flag[horizontal] > 0) {
                dqdt[index] = humidity_factor[horizontal] *
                    local_dq[index] / deltsec2 - dqt[index];
                dtdt[index] = temperature_factor[horizontal] *
                    local_dt[index] / deltsec2 + dtdt[index];
                layer_convective_rain[index] =
                    temperature_factor[horizontal] * local_dt[index] *
                    sigma_thickness[level] / latent_over_cp[index] /
                    deltsec2;
            }
            if (layer_convective_rain[index] > 0.0f) {
                convective_layer_flag[index] = 1;
                convective_layer_count[horizontal] += 1;
            }
            if (sigma[level] < 0.4f &&
                convective_layer_flag[index] == 1 &&
                top_level[horizontal] == level + 1) {
                convective_layer_flag[index] = 2;
            }
        }
    }

    if (shallow_convection == 1) {
        status = rain_shallow_convection(
            nhor, nlev, mix_momentum, dissipative_heating,
            shallow_diffusion, shallow_top_pressure, gravity, deltsec2,
            gascon, kap, rdbrv, ra1, ra2, tmelt, ra4, cpd, cpv_cpd_minus1, sigma,
            sigma_half, sigma_thickness, surface_pressure, du, dv, dq,
            dudt, dvdt, environment_t, environment_q, shallow_flag,
            lift_level, top_level, shallow_dt, shallow_dq);
        if (status != 0) {
            free(integer_storage);
            free(column_storage);
            free(field_storage);
            return status;
        }
        for (level = 0; level < nlev; ++level) {
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                dtdt[index] += shallow_dt[index];
                dqdt[index] += shallow_dq[index];
            }
        }
    }

    if (mix_momentum == 1) {
        for (level = 0; level < nlev; ++level) {
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                environment_u[index] = du[index] + dudt[index] * deltsec2;
                environment_v[index] = dv[index] + dvdt[index] * deltsec2;
            }
        }
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            momentum_factor[horizontal] = 0.0f;
            if (temperature_energy[horizontal] +
                    humidity_energy[horizontal] > 0.0f) {
                momentum_factor[horizontal] = rain_min(
                    1.0f, available_energy[horizontal] /
                        (temperature_energy[horizontal] +
                         humidity_energy[horizontal]));
            }
        }
        for (level = 0; level < nlev; ++level) {
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                if (lift_level[horizontal] >= level + 1 &&
                    top_level[horizontal] <= level + 1 &&
                    available_energy[horizontal] > 0.0f &&
                    shallow_flag[horizontal] == 0) {
                    parcel_u[horizontal] +=
                        environment_u[index] * sigma_thickness[level];
                    parcel_v[horizontal] +=
                        environment_v[index] * sigma_thickness[level];
                    mixed_sigma[horizontal] += sigma_thickness[level];
                }
            }
        }
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            if (mixed_sigma[horizontal] > 0.0f) {
                parcel_u[horizontal] /= mixed_sigma[horizontal];
                parcel_v[horizontal] /= mixed_sigma[horizontal];
            }
        }
        for (level = 0; level < nlev; ++level) {
            for (horizontal = 0; horizontal < nhor; ++horizontal) {
                const size_t index =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;
                if (lift_level[horizontal] >= level + 1 &&
                    top_level[horizontal] <= level + 1 &&
                    available_energy[horizontal] > 0.0f &&
                    shallow_flag[horizontal] == 0) {
                    local_du[index] = momentum_factor[horizontal] *
                        (parcel_u[horizontal] - environment_u[index]);
                    local_dv[index] = momentum_factor[horizontal] *
                        (parcel_v[horizontal] - environment_v[index]);
                    dudt[index] += local_du[index] / deltsec2;
                    dvdt[index] += local_dv[index] / deltsec2;
                }
            }
        }
    }

    free(integer_storage);
    free(column_storage);
    free(field_storage);
    return 0;
}

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
    float *surface_snow)
{
    const size_t field_count = (size_t)nhor * (size_t)nlev;
    float *storage = (float *)calloc(
        (size_t)nhor * 4U + field_count * 3U, sizeof(float));
    float *convective_rain;
    float *large_scale_rain;
    float *convective_snow;
    float *large_scale_snow;
    float *temperature_tendency;
    float *rain_humidity_tendency;
    float *snow_humidity_tendency;
    int32_t level;
    int32_t horizontal;

    if (storage == NULL) {
        return 1;
    }
    convective_rain = storage;
    large_scale_rain = convective_rain + nhor;
    convective_snow = large_scale_rain + nhor;
    large_scale_snow = convective_snow + nhor;
    temperature_tendency = large_scale_snow + nhor;
    rain_humidity_tendency = temperature_tendency + field_count;
    snow_humidity_tendency = rain_humidity_tendency + field_count;

    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const float layer_cp =
                (ls - lv) /
                (cpd * (1.0f + cpv_cpd_minus1 * dq[index]));
            convective_rain[horizontal] =
                convective_rain[horizontal] +
                layer_convective_rain[index] -
                layer_convective_snow[index];
            large_scale_rain[horizontal] =
                large_scale_rain[horizontal] +
                layer_large_scale_rain[index] -
                layer_large_scale_snow[index];
            convective_snow[horizontal] += layer_convective_snow[index];
            large_scale_snow[horizontal] += layer_large_scale_snow[index];

            if (convective_snow[horizontal] +
                    large_scale_snow[horizontal] > 0.0f &&
                dt[index] > tmelt) {
                temperature_tendency[index] =
                    -(convective_snow[horizontal] +
                      large_scale_snow[horizontal]) /
                    sigma_thickness[level] * layer_cp;
                convective_rain[horizontal] =
                    convective_snow[horizontal] +
                    convective_rain[horizontal];
                large_scale_rain[horizontal] =
                    large_scale_snow[horizontal] +
                    large_scale_rain[horizontal];
                convective_snow[horizontal] = 0.0f;
                large_scale_snow[horizontal] = 0.0f;
            } else if (convective_rain[horizontal] +
                           large_scale_rain[horizontal] > 0.0f &&
                       dt[index] <= tmelt) {
                temperature_tendency[index] =
                    (convective_rain[horizontal] +
                     large_scale_rain[horizontal]) /
                    sigma_thickness[level] * layer_cp;
                convective_snow[horizontal] =
                    convective_rain[horizontal] +
                    convective_snow[horizontal];
                large_scale_snow[horizontal] =
                    large_scale_rain[horizontal] +
                    large_scale_snow[horizontal];
                convective_rain[horizontal] = 0.0f;
                large_scale_rain[horizontal] = 0.0f;
            }

            if (level == nlev - 1) {
                const size_t surface =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)nlev;
                if (convective_snow[horizontal] +
                        large_scale_snow[horizontal] > 0.0f &&
                    dt[surface] > tmelt) {
                    temperature_tendency[index] =
                        -(convective_snow[horizontal] +
                          large_scale_snow[horizontal]) /
                        sigma_thickness[level] * layer_cp +
                        temperature_tendency[index];
                    convective_rain[horizontal] =
                        convective_snow[horizontal] +
                        convective_rain[horizontal];
                    large_scale_rain[horizontal] =
                        large_scale_snow[horizontal] +
                        large_scale_rain[horizontal];
                    convective_snow[horizontal] = 0.0f;
                    large_scale_snow[horizontal] = 0.0f;
                } else if (convective_rain[horizontal] +
                               large_scale_rain[horizontal] > 0.0f &&
                           dt[surface] <= tmelt) {
                    temperature_tendency[index] =
                        (convective_rain[horizontal] +
                         large_scale_rain[horizontal]) /
                        sigma_thickness[level] * layer_cp +
                        temperature_tendency[index];
                    convective_snow[horizontal] =
                        convective_rain[horizontal] +
                        convective_snow[horizontal];
                    large_scale_snow[horizontal] =
                        large_scale_rain[horizontal] +
                        large_scale_snow[horizontal];
                    convective_rain[horizontal] = 0.0f;
                    large_scale_rain[horizontal] = 0.0f;
                }
            }

            if (evaporate_precipitation == 1) {
                const float predicted_q =
                    dq[index] + dqdt[index] * deltsec2;
                const float predicted_t =
                    dt[index] + dtdt[index] * deltsec2;
                float saturation = rain_saturation_mixing_ratio(
                    predicted_t, surface_pressure[horizontal], sigma[level],
                    rdbrv, ra1, ra2, tmelt, ra4);
                const float correction =
                    rain_mixing_ratio_correction(saturation, rdbrv);
                float latent_over_cp;
                float evaporated;
                float denominator;

                saturation = saturation * correction;
                if (saturation < 0.0f) {
                    saturation = 0.0f;
                }
                denominator = predicted_t - ra4;
                denominator = denominator * denominator;

#define RAIN_EVAPORATE(amount, latent, tendency)                         \
                do {                                                     \
                    if (predicted_q < saturation && (amount) > 0.0f) {   \
                        latent_over_cp = (latent) /                       \
                            (cpd * (1.0f + cpv_cpd_minus1 * dq[index]));            \
                        evaporated = gamma *                              \
                            (saturation - predicted_q) *                  \
                            sigma_thickness[level] / deltsec2 /           \
                            (1.0f + latent_over_cp * ra2 *                \
                             (tmelt - ra4) * saturation * correction /   \
                             denominator);                               \
                        evaporated = rain_min(evaporated, (amount));      \
                        (amount) = (amount) - evaporated;                 \
                        evaporated =                                     \
                            evaporated / sigma_thickness[level];          \
                        temperature_tendency[index] =                     \
                            temperature_tendency[index] -                 \
                            evaporated * latent_over_cp;                  \
                        (tendency)[index] =                               \
                            (tendency)[index] + evaporated;               \
                    }                                                     \
                } while (0)

                RAIN_EVAPORATE(
                    large_scale_rain[horizontal], lv,
                    rain_humidity_tendency);
                RAIN_EVAPORATE(
                    convective_rain[horizontal], lv,
                    rain_humidity_tendency);
                RAIN_EVAPORATE(
                    large_scale_snow[horizontal], ls,
                    snow_humidity_tendency);
                RAIN_EVAPORATE(
                    convective_snow[horizontal], ls,
                    snow_humidity_tendency);
#undef RAIN_EVAPORATE
            }
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal) {
        surface_convective_rain[horizontal] =
            (convective_rain[horizontal] +
             convective_snow[horizontal]) *
            surface_pressure[horizontal] / gravity / 1000.0f;
        surface_large_scale_rain[horizontal] =
            (large_scale_rain[horizontal] +
             large_scale_snow[horizontal]) *
            surface_pressure[horizontal] / gravity / 1000.0f;
        surface_snow[horizontal] =
            (convective_snow[horizontal] +
             large_scale_snow[horizontal]) *
            surface_pressure[horizontal] / gravity / 1000.0f;
    }
    for (level = 0; level < nlev; ++level) {
        for (horizontal = 0; horizontal < nhor; ++horizontal) {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            dqdt[index] = dqdt[index] +
                rain_humidity_tendency[index] +
                snow_humidity_tendency[index];
            dtdt[index] = dtdt[index] + temperature_tendency[index];
        }
    }

    free(storage);
    return 0;
}

#undef RAIN_AT
