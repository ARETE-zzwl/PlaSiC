/*
 * 在线输出的运行时实现：写出进度记录与二进制帧文件。
 * Runtime implementation of live output: progress records and binary frame
 * files.
 *
 * 两类产物 / Two output kinds:
 *   progress_output：progress.txt 进度文件，每次写入追加一行 JSON，只包含
 *                    completed_steps、total_steps、date、elapsed_seconds 和
 *                    steps_per_second 五个字段。
 *                    progress.txt progress file; each write appends one JSON
 *                    line holding only completed_steps, total_steps, date,
 *                    elapsed_seconds and steps_per_second.
 *   frame_output：   二维/三维格点场的二进制帧，按 frame_interval 降采样写出。
 *                    Binary stream of 2-D/3-D grid-point fields, written every
 *                    frame_interval steps.
 *
 * 并行语义 / Parallel semantics: 每个 rank 累加自己的本地诊断并提供纬度条带；
 * mp_allgather_grid() 拼出全球场后，只有根进程持有文件句柄并写盘，写盘结果再广播
 * 给所有 rank，保证大家一致地继续或停止。
 * Every rank accumulates its local diagnostics and contributes its latitude
 * band to mp_allgather_grid(); only the root holds file handles and writes
 * files, and the result is broadcast so that all ranks always agree on
 * whether to continue or stop.
 */

#include "runtime_internal.h"

/*
 * 检查全部预报数组是否都是有限值（无 NaN/Inf）。
 * Check that every prognostic array holds only finite values (no NaN/Inf).
 */
int runtime_state_is_finite(void)
{
    size_t index;

    for (index = 0; index < sizeof(plasic_st) / sizeof(plasic_st[0]);
         ++index)
    {
        if (!isfinite(plasic_st[index]) ||
            !isfinite(plasic_sd[index]) ||
            !isfinite(plasic_sz[index]) ||
            !isfinite(plasic_sq[index]))
        {
            return 0;
        }
    }
    for (index = 0; index < sizeof(plasic_sp) / sizeof(plasic_sp[0]);
         ++index)
    {
        if (!isfinite(plasic_sp[index]))
        {
            return 0;
        }
    }
    return 1;
}

static void reset_output_accumulators(void)
{
    size_t horizontal;

    memset(plasic_aprl, 0, sizeof(plasic_aprl));
    memset(plasic_aprc, 0, sizeof(plasic_aprc));
    memset(plasic_aprs, 0, sizeof(plasic_aprs));
    memset(plasic_aevap, 0, sizeof(plasic_aevap));
    memset(plasic_ashfl, 0, sizeof(plasic_ashfl));
    memset(plasic_alhfl, 0, sizeof(plasic_alhfl));
    memset(plasic_asmelt, 0, sizeof(plasic_asmelt));
    memset(plasic_asndch, 0, sizeof(plasic_asndch));
    memset(plasic_acc, 0, sizeof(plasic_acc));
    memset(plasic_assol, 0, sizeof(plasic_assol));
    memset(plasic_asthr, 0, sizeof(plasic_asthr));
    memset(plasic_atsol, 0, sizeof(plasic_atsol));
    memset(plasic_atthr, 0, sizeof(plasic_atthr));
    memset(plasic_ataux, 0, sizeof(plasic_ataux));
    memset(plasic_atauy, 0, sizeof(plasic_atauy));
    memset(plasic_atsolu, 0, sizeof(plasic_atsolu));
    memset(plasic_assolu, 0, sizeof(plasic_assolu));
    memset(plasic_asthru, 0, sizeof(plasic_asthru));
    memset(plasic_aqvi, 0, sizeof(plasic_aqvi));
    memset(plasic_atsa, 0, sizeof(plasic_atsa));
    memset(plasic_ats0, 0, sizeof(plasic_ats0));
    memset(plasic_atsama, 0, sizeof(plasic_atsama));
    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        plasic_atsami[horizontal] = 1.0e10f;
    }
    plasic_naccuout = 0;
}

void runtime_accumulate_output_fields(void)
{
    const size_t surface_offset =
        (size_t)PLASIC_NHOR * (size_t)PLASIC_NLEV;
    size_t horizontal;

    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        plasic_aprl[horizontal] += plasic_dprl[horizontal];
        plasic_aprc[horizontal] += plasic_dprc[horizontal];
        plasic_aprs[horizontal] += plasic_dprs[horizontal];
        plasic_aevap[horizontal] += plasic_devap[horizontal];
        plasic_ashfl[horizontal] += plasic_dshfl[horizontal];
        plasic_alhfl[horizontal] += plasic_dlhfl[horizontal];
        plasic_asmelt[horizontal] += plasic_dsmelt[horizontal];
        plasic_asndch[horizontal] += plasic_dsndch[horizontal];
        plasic_acc[horizontal] +=
            plasic_dcc[surface_offset + horizontal];
        plasic_assol[horizontal] +=
            plasic_dswfl[surface_offset + horizontal];
        plasic_asthr[horizontal] +=
            plasic_dlwfl[surface_offset + horizontal];
        plasic_atsol[horizontal] += plasic_dswfl[horizontal];
        plasic_atthr[horizontal] += plasic_dlwfl[horizontal];
        plasic_assolu[horizontal] +=
            plasic_dfu[surface_offset + horizontal];
        plasic_asthru[horizontal] +=
            plasic_dftu[surface_offset + horizontal];
        plasic_atsolu[horizontal] += plasic_dfu[horizontal];
        plasic_ataux[horizontal] += plasic_dtaux[horizontal];
        plasic_atauy[horizontal] += plasic_dtauy[horizontal];
        plasic_aqvi[horizontal] += plasic_dqvi[horizontal];
        plasic_atsa[horizontal] += plasic_dtsa[horizontal];
        plasic_ats0[horizontal] +=
            plasic_dt[surface_offset + horizontal];
        plasic_atsami[horizontal] =
            fminf(plasic_atsami[horizontal],
                  plasic_dtsa[horizontal]);
        plasic_atsama[horizontal] =
            fmaxf(plasic_atsama[horizontal],
                  plasic_dtsa[horizontal]);
    }
    ++plasic_naccuout;
    if (plasic_nafter > 0 &&
        plasic_nstep % plasic_nafter == 0)
    {
        reset_output_accumulators();
    }
}

/**
 * 获取当前墙钟时间（秒）。
 * 仅用于统计实际运行耗时，不参与模拟时间推进。
 * return 当前 UTC 时间对应的秒数；获取失败返回 0.0。
 *
 * Return the current wall-clock time in seconds.
 * It is only used to measure real elapsed time and never advances model time.
 * Returns the UTC time in seconds, or 0.0 when the time source fails.
 */
static double wall_time_seconds(void)
{
    struct timespec now;

    if (timespec_get(&now, TIME_UTC) != TIME_UTC)
    {
        return 0.0;
    }
    return (double)now.tv_sec + (double)now.tv_nsec * 1.0e-9;
    // tv_sec = time value seconds
    // tv_nsec = time value nanoseconds
}

static int write_progress_record(
    const plasic_runtime_options *options, int32_t completed_steps)
{
    double elapsed;
    double rate; /* 平均计算速度，单位为 step/s。 */
                 /* Mean computation rate in steps per second. */

    if (!progress_enabled || !mp_is_root())
    {
        return PLASIC_RUNTIME_OK;
    }

    elapsed = wall_time_seconds() - integration_start_time;

    rate = elapsed > 0.0 ? (double)completed_steps / elapsed : 0.0;

    fprintf(
        progress_stream,
        "{\"completed_steps\":%d,\"total_steps\":%d,"
        "\"date\":{\"year\":%d,\"month\":%d,\"day\":%d,"
        "\"hour\":%d,\"minute\":%d},"
        "\"elapsed_seconds\":%.9g,\"steps_per_second\":%.9g}\n",
        completed_steps, options->steps,
        plasic_ndatim[0], plasic_ndatim[1],
        plasic_ndatim[2], plasic_ndatim[3],
        plasic_ndatim[4], elapsed, rate);

    return fflush(progress_stream) == 0 ? PLASIC_RUNTIME_OK : PLASIC_RUNTIME_IO_FAILED;
}

static int write_progress_step(
    const plasic_runtime_options *options, int32_t completed_steps)
{
    /*
     * 每完成一个完整模式时间步检查一次是否追加进度记录。只有到达
     * progress_interval 的整数倍或完成本次运行最后一步时才真正写文件
     *
     * After each complete model time step, decide whether to append a progress
     * record. The file is written only on multiples of progress_interval or on
     * the final step of this run
     */

    if (!progress_enabled ||
        (completed_steps % options->progress_interval != 0 &&
         completed_steps != options->steps))
    {
        return PLASIC_RUNTIME_OK;
    }

    return write_progress_record(options, completed_steps);
}

static int write_frame_field(
    int32_t variable_id, const char *name, const char *unit,
    const float *local_grid_values, int32_t layer_count)
{
    /*
     * 将一个格点变量写成 App 可读取的二进制帧。
     *
     * MPI 将全球高斯网格按纬度条带分给各 rank，而帧文件要求每条记录包含完整
     * 全球场。因此所有 rank 先把各自的局部格点条带拼成 global_grid_values，
     * 然后只由根进程写盘，最后把写入状态广播给其余进程。
     *
     * 参数含义：
     *   variable_id：plasic_stream_variable 中定义的标准变量编号；
     *   name/unit：随帧写入的变量名和单位元数据；
     *   local_grid_values：本 rank 的格点数组，布局为 [layer][local point]；
     *   layer_count：此次 payload 包含的垂直层数，二维地表场传 1。
     *
     * 本函数只整理和输出数据，不修改模式预报状态。
     *
     * Write one grid-point variable as a binary frame the App can read.
     *
     * MPI splits the global Gaussian grid into latitude bands, one per rank,
     * while each frame record must hold the complete global field. Every rank
     * therefore contributes its local band to global_grid_values first, then
     * only the root writes the file, and finally the write status is
     * broadcast back to the other ranks.
     *
     * Parameters:
     *   variable_id: standard variable id declared in plasic_stream_variable;
     *   name/unit: variable name and unit metadata written with the frame;
     *   local_grid_values: this rank's grid array, laid out as
     *                      [layer][local point];
     *   layer_count: vertical layers in this payload, 1 for 2-D surface
     *                fields.
     *
     * This function only reorganizes and emits data; it never modifies the
     * prognostic state.
     */
    float *global_grid_values; /* 汇集后的 [layer][global point] 完整全球场。 */
    int status;

    if (local_grid_values == NULL || layer_count <= 0)
    {
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    global_grid_values = (float *)malloc(
        (size_t)PLASIC_NUGP * (size_t)layer_count *
        sizeof(*global_grid_values));
    if (global_grid_values == NULL)
    {
        return PLASIC_RUNTIME_IO_FAILED;
    }

    /*
     * MPI 模式下收集并按全局纬度顺序拼接所有局部条带；串行模式下等价于
     * 把 local_grid_values 复制到 global_grid_values。
     *
     * Under MPI the local bands are collected and concatenated in global
     * latitude order; in serial mode this is equivalent to copying
     * local_grid_values into global_grid_values.
     */
    status = mp_allgather_grid(
                 global_grid_values, local_grid_values,
                 PLASIC_NUGP, PLASIC_NHOR, layer_count) == 0
                 ? PLASIC_RUNTIME_OK
                 : PLASIC_RUNTIME_NUMERICAL_FAILURE;

    if (status == PLASIC_RUNTIME_OK && mp_is_root())
    {
        status = plasic_stream_write_frame(
                     &frame_writer, variable_id, name, unit,
                     (int64_t)plasic_nstep, plasic_ndatim,
                     layer_count == 1 ? 0 : -1, layer_count,
                     global_grid_values) ==
                         PLASIC_STREAM_OK
                     ? PLASIC_RUNTIME_OK
                     : PLASIC_RUNTIME_IO_FAILED;
    }

    status = runtime_broadcast_status(status);

    free(global_grid_values);
    return status;
}

static int write_live_frames(void)
{
    /*
     * 把内部变量及各物理圈层的诊断投影为 App 可读的标准字段。这里是“输出适配层”，
     * 不应改变积分状态；派生量（风速、相对湿度、位势高度等）全部写到工作区。
     *
     * Project the internal variables and diagnostics of each physical
     * component into standard fields the App can read. This is an output
     * adapter layer and must not change the integration state; derived
     * quantities (wind speed, relative humidity, geopotential height, ...)
     * are all written into the work arrays.
     */
    const float precipitation_to_mm_per_hour = 3.6e6f;
    const float *physics_values;
    int32_t layers;
    int32_t level;
    int status;
    size_t horizontal;

#define WRITE_FRAME(id, name, unit, values, count)    \
    do                                                \
    {                                                 \
        status = write_frame_field(                   \
            (id), (name), (unit), (values), (count)); \
        if (status != PLASIC_RUNTIME_OK)              \
        {                                             \
            return status;                            \
        }                                             \
    } while (0)

    WRITE_FRAME(
        PLASIC_STREAM_AIR_TEMPERATURE, "air_temperature", "K",
        plasic_dt, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_SPECIFIC_HUMIDITY, "specific_humidity",
        "kg kg-1", plasic_dq, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_CLOUD_COVER, "cloud_cover", "1",
        plasic_dcc, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_ZONAL_WIND, "zonal_wind", "m s-1",
        plasic_du, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_MERIDIONAL_WIND, "meridional_wind", "m s-1",
        plasic_dv, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_VERTICAL_VELOCITY, "pressure_velocity", "Pa s-1",
        plasic_dw, PLASIC_NLEV);
    /* gd/gz 分别为无量纲的 D/Ω、绝对涡度/Ω 场；下面恢复量纲并拆出行星涡度。 */
    /* gd/gz are nondimensional D/Omega and absolute-zeta/Omega fields; the loop below restores dimensions and removes the planetary vorticity. */
    for (level = 0; level < PLASIC_NLEV; ++level)
    {
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            const size_t cell = horizontal + (size_t)level * PLASIC_NHOR;
            const size_t latitude = horizontal / PLASIC_NLON;
            const float planetary_vorticity =
                plasic_plavor * sqrtf(1.5f) * (float)plasic_sid[latitude];

            runtime_fields->physical_divergence[cell] =
                plasic_ww * plasic_gd[cell];
            runtime_fields->relative_vorticity[cell] =
                plasic_ww * (plasic_gz[cell] - planetary_vorticity);
        }
    }
    WRITE_FRAME(
        PLASIC_STREAM_VORTICITY, "relative_vorticity", "s-1",
        runtime_fields->relative_vorticity, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_DIVERGENCE, "divergence", "s-1",
        runtime_fields->physical_divergence, PLASIC_NLEV);

    /* 风速为水平风矢量的模；派生量写入工作区，不改动预报风场。 */
    /* Wind speed is the magnitude of the horizontal wind vector; the derived values go to the work array and never modify the prognostic winds. */
    for (horizontal = 0;
         horizontal < (size_t)PLASIC_NHOR * (size_t)PLASIC_NLEV;
         ++horizontal)
    {
        runtime_fields->wind_speed[horizontal] = sqrtf(
            plasic_du[horizontal] * plasic_du[horizontal] +
            plasic_dv[horizontal] * plasic_dv[horizontal]);
    }
    WRITE_FRAME(
        PLASIC_STREAM_WIND_SPEED, "wind_speed", "m s-1",
        runtime_fields->wind_speed, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_SURFACE_TEMPERATURE, "surface_temperature", "K",
        plasic_dt + (size_t)PLASIC_NHOR * (size_t)PLASIC_NLEV,
        1);

    /* 降水诊断单位为 m/s 等效水柱，乘以 3.6e6（即 1000 * 3600）换算为 mm/hr。 */
    /*
     * Precipitation diagnostics are in m/s water equivalent and are converted
     * to mm/hr by the 3.6e6 factor (1000 * 3600).
     */
    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        runtime_fields->total_precipitation_rate[horizontal] =
            (plasic_dprl[horizontal] +
             plasic_dprc[horizontal]) *
            precipitation_to_mm_per_hour;
        runtime_fields->convective_precipitation_rate[horizontal] =
            plasic_dprc[horizontal] * precipitation_to_mm_per_hour;
        runtime_fields->large_scale_precipitation_rate[horizontal] =
            plasic_dprl[horizontal] * precipitation_to_mm_per_hour;
        runtime_fields->snowfall_rate[horizontal] =
            plasic_dprs[horizontal] * precipitation_to_mm_per_hour;
    }
    WRITE_FRAME(
        PLASIC_STREAM_PRECIPITATION, "precipitation", "mm/hr",
        runtime_fields->total_precipitation_rate, 1);

    WRITE_FRAME(
        PLASIC_STREAM_SURFACE_PRESSURE, "surface_pressure", "Pa",
        plasic_dp, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SURFACE_AIR_TEMPERATURE, "surface_air_temperature",
        "K", plasic_dtsa, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SURFACE_HUMIDITY, "surface_relative_humidity", "1",
        plasic_dwetfac, 1);
    WRITE_FRAME(
        PLASIC_STREAM_EVAPORATION, "evaporation", "m s-1",
        plasic_devap, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SENSIBLE_HEAT_FLUX, "sensible_heat_flux", "W m-2",
        plasic_dshfl, 1);
    WRITE_FRAME(
        PLASIC_STREAM_LATENT_HEAT_FLUX, "latent_heat_flux", "W m-2",
        plasic_dlhfl, 1);
    WRITE_FRAME(
        PLASIC_STREAM_ALBEDO, "albedo", "1", plasic_dalb, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SNOW_DEPTH, "snow_water_equivalent", "m",
        plasic_dsnow, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SOIL_TEMPERATURE, "soil_temperature", "K",
        plasic_dtsoil, 1);
    WRITE_FRAME(
        PLASIC_STREAM_CONVECTIVE_PRECIPITATION,
        "convective_precipitation", "mm/hr",
        runtime_fields->convective_precipitation_rate, 1);
    WRITE_FRAME(
        PLASIC_STREAM_LARGE_SCALE_PRECIPITATION,
        "large_scale_precipitation", "mm/hr",
        runtime_fields->large_scale_precipitation_rate, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SNOWFALL, "snowfall", "mm/hr",
        runtime_fields->snowfall_rate, 1);
    WRITE_FRAME(
        PLASIC_STREAM_ZONAL_WIND_STRESS, "zonal_wind_stress", "N m-2",
        plasic_dtaux, 1);
    WRITE_FRAME(
        PLASIC_STREAM_MERIDIONAL_WIND_STRESS, "meridional_wind_stress",
        "N m-2", plasic_dtauy, 1);
    WRITE_FRAME(
        PLASIC_STREAM_LAND_MASK, "land_mask", "1",
        plasic_dls, 1);
    WRITE_FRAME(
        PLASIC_STREAM_MIXED_LAYER_DEPTH, "mixed_layer_depth", "m",
        plasic_dmld, 1);
    WRITE_FRAME(
        PLASIC_STREAM_LIQUID_WATER, "liquid_water", "kg kg-1",
        plasic_dql, PLASIC_NLEV);

    /* 总云量采用最大—随机重叠：1 - Π_k (1 - cc_k)。 */
    /* Total cloud cover uses maximum-random overlap: 1 - Π_k (1 - cc_k). */
    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        float clear_sky = 1.0f;

        for (level = 0; level < PLASIC_NLEV; ++level)
        {
            clear_sky *= 1.0f - plasic_dcc[horizontal + (size_t)level * PLASIC_NHOR];
        }
        runtime_fields->total_cloud_cover[horizontal] = 1.0f - clear_sky;
    }
    WRITE_FRAME(
        PLASIC_STREAM_TOTAL_CLOUD_COVER, "total_cloud_cover", "1",
        runtime_fields->total_cloud_cover, 1);

    for (level = 0; level < PLASIC_NLEV; ++level)
    {
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            const size_t cell =
                horizontal + (size_t)level * PLASIC_NHOR;
            const float air_temperature = plasic_dt[cell];
            const float level_pressure =
                plasic_dp[horizontal] * plasic_sigma[level];
            const float saturation_mixing_ratio = fminf(
                plasic_rdbrv * plasic_ra1 *
                    expf(plasic_ra2 *
                         (air_temperature - plasic_tmelt) /
                         (air_temperature - plasic_ra4)) /
                    level_pressure,
                plasic_rdbrv);
            const float saturation_specific_humidity =
                saturation_mixing_ratio /
                (1.0f - (1.0f / plasic_rdbrv - 1.0f) *
                            saturation_mixing_ratio);

            runtime_fields->relative_humidity[cell] =
                plasic_dq[cell] / saturation_specific_humidity;
        }
    }
    WRITE_FRAME(
        PLASIC_STREAM_RELATIVE_HUMIDITY, "relative_humidity", "1",
        runtime_fields->relative_humidity, PLASIC_NLEV);

    /*
     * 位势高度由静力关系从地面 σ=1 向上积分，并加上初始化得到的地形位势。
     * sigmah[level] 是该层下边界，sigma[level] 是层中心。
     *
     * Geopotential height is integrated upward from the surface (sigma = 1)
     * with the hydrostatic relation, then the surface geopotential from
     * initialization is added. sigmah[level] is the lower boundary of the
     * layer and sigma[level] is the layer centre.
     */
    for (level = 0; level < PLASIC_NLEV; ++level)
    {
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            float hydrostatic_temperature_integral = 0.0f;
            int32_t lower_level;

            for (lower_level = level + 1;
                 lower_level < PLASIC_NLEV; ++lower_level)
            {
                hydrostatic_temperature_integral += plasic_ct *
                                                    (plasic_gt[horizontal + (size_t)lower_level * PLASIC_NHOR] +
                                                     runtime_fields->nondimensional_reference_temperature[lower_level]) *
                                                    logf(plasic_sigmah[lower_level] /
                                                         plasic_sigmah[lower_level - 1]);
            }
            hydrostatic_temperature_integral += plasic_ct *
                                                (plasic_gt[horizontal + (size_t)level * PLASIC_NHOR] +
                                                 runtime_fields->nondimensional_reference_temperature[level]) *
                                                logf(plasic_sigmah[level] /
                                                     plasic_sigma[level]);
            runtime_fields->geopotential_height[horizontal + (size_t)level * PLASIC_NHOR] =
                (runtime_fields->surface_geopotential[horizontal] +
                 plasic_gascon * hydrostatic_temperature_integral) /
                plasic_ga;
        }
    }
    WRITE_FRAME(
        PLASIC_STREAM_GEOPOTENTIAL_HEIGHT, "geopotential_height", "m",
        runtime_fields->geopotential_height, PLASIC_NLEV);
    WRITE_FRAME(
        PLASIC_STREAM_SEA_ICE_FRACTION, "sea_ice_fraction", "1",
        plasic_dicec, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SEA_ICE_THICKNESS, "sea_ice_thickness", "m",
        plasic_diced, 1);
    WRITE_FRAME(
        PLASIC_STREAM_SOIL_WATER, "soil_water", "m",
        plasic_dwatc, 1);

    status = plasic_physics_frame_field(
        PLASIC_STREAM_SEA_SURFACE_TEMPERATURE, &physics_values, &layers);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    WRITE_FRAME(
        PLASIC_STREAM_SEA_SURFACE_TEMPERATURE,
        "sea_surface_temperature", "K", physics_values, layers);

    /*
     * 三维深海温度与盐度为可选输出：深海模块未启用或尚未初始化时，
     * plasic_physics_frame_field 返回 BAD_ARGUMENT，此时跳过而不写帧。
     *
     * Three-dimensional deep-ocean temperature and salinity are optional
     * outputs: when the deep-ocean module is disabled or not initialized yet,
     * plasic_physics_frame_field returns BAD_ARGUMENT and the frames are
     * skipped instead of being written.
     */
    status = plasic_physics_frame_field(
        PLASIC_STREAM_DEEP_OCEAN_TEMPERATURE, &physics_values, &layers);
    if (status == PLASIC_RUNTIME_OK)
    {
        WRITE_FRAME(
            PLASIC_STREAM_DEEP_OCEAN_TEMPERATURE,
            "deep_ocean_temperature", "K", physics_values, layers);
    }
    else if (status != PLASIC_RUNTIME_BAD_ARGUMENT)
    {
        return status;
    }

    status = plasic_physics_frame_field(
        PLASIC_STREAM_DEEP_OCEAN_SALINITY, &physics_values, &layers);
    if (status == PLASIC_RUNTIME_OK)
    {
        WRITE_FRAME(
            PLASIC_STREAM_DEEP_OCEAN_SALINITY,
            "deep_ocean_salinity", "g kg-1", physics_values, layers);
    }
    else if (status != PLASIC_RUNTIME_BAD_ARGUMENT)
    {
        return status;
    }

#undef WRITE_FRAME
    return PLASIC_RUNTIME_OK;
}

int runtime_write_step_outputs(
    const plasic_runtime_options *options, int32_t completed_steps)
{
    /*
     * 写出一个完整模式时间步结束后的在线输出。
     *
     * 输出分两类：
     *   1. 供 GUI/监控程序读取的 progress.txt 进度记录；
     *   2. 较大的多变量二进制帧，按 frame_interval 降采样写出。
     *
     * Write the live outputs at the end of one complete model time step.
     *
     * There are two output kinds:
     *   1. progress.txt records read by the GUI/monitor;
     *   2. larger multi-variable binary frames, written at the frame_interval
     *      downsampling rate.
     */
    int status;

    status = write_progress_step(options, completed_steps);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    if (frames_enabled &&
        (completed_steps % options->frame_interval == 0 ||
         completed_steps == options->steps))
    {
        status = write_live_frames();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
    }

    return PLASIC_RUNTIME_OK;
}

/*
 * 帧文件头中声明的海洋垂直层数：三维深海海洋编译进模型时为其实
 * 际层数，否则为单层混合层海洋的 1 层。App 等下游读取方据此判断
 * 流中是否可能包含三维海洋变量；帧本身仍逐帧携带自己的 layers。
 *
 * Ocean level count advertised in the frame-stream file header: the
 * actual deep-ocean layer count when the three-dimensional ocean is
 * compiled in, otherwise 1 for the single mixed layer. Downstream
 * readers use it to tell whether 3-D ocean variables can appear in
 * the stream; each frame still carries its own layer count.
 */
static int32_t stream_ocean_levels(void)
{
#if defined(PLASIC_DEEP_OCEAN) && PLASIC_DEEP_OCEAN && \
    defined(PLASIC_DEEP_NLEV)
    return PLASIC_DEEP_NLEV;
#else
    return 1;
#endif
}

int runtime_open_live_outputs(const plasic_runtime_options *options)
{
    int status = PLASIC_RUNTIME_OK;

    progress_enabled = options->progress_output != NULL;
    frames_enabled = options->frame_output != NULL;

    if (mp_is_root() && progress_enabled)
    {
        progress_stream = fopen(options->progress_output, "w");
        if (progress_stream == NULL)
        {
            runtime_set_error(
                "cannot create progress file '%s'",
                options->progress_output);
            status = PLASIC_RUNTIME_IO_FAILED;
        }
    }

    status = runtime_broadcast_status(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    if (mp_is_root() && frames_enabled)
    {
        if (plasic_stream_open(
                &frame_writer, options->frame_output, PLASIC_NLON,
                PLASIC_NLAT, PLASIC_NLEV, stream_ocean_levels()) !=
            PLASIC_STREAM_OK)
        {
            runtime_set_error(
                "cannot create frame stream '%s'",
                options->frame_output);
            status = PLASIC_RUNTIME_IO_FAILED;
        }
    }

    status = runtime_broadcast_status(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    integration_start_time = wall_time_seconds();

    return PLASIC_RUNTIME_OK;
}

int runtime_close_step_outputs(const plasic_runtime_options *options)
{
    int status = PLASIC_RUNTIME_OK;

    /* 参数仅供统一的运行时输出接口使用，本函数不再读取其中字段。 */
    /* The parameter is kept for the common runtime output interface; this
     * function no longer reads any of its fields. */
    (void)options;

    if (mp_is_root() && progress_stream != NULL)
    {
        if (fclose(progress_stream) != 0)
        {
            status = PLASIC_RUNTIME_IO_FAILED;
        }
        progress_stream = NULL;
    }

    if (mp_is_root() && frame_writer.stream != NULL)
    {
        if (plasic_stream_close(&frame_writer) != PLASIC_STREAM_OK)
        {
            status = PLASIC_RUNTIME_IO_FAILED;
        }
    }

    return runtime_broadcast_status(status);
}
