#include "runtime_internal.h"

static int synthesize_adiabatic_fields(void)
{
    /*
     * Spectral → grid synthesis
     */
    int status = sht_vortdiv_to_wind(
        PLASIC_NLEV, plasic_plavor, plasic_sd, plasic_sz,
        plasic_gu, plasic_gv);

    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(
            PLASIC_NLEV, plasic_sd, plasic_gd);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(
            PLASIC_NLEV, plasic_st, plasic_gt);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(
            PLASIC_NLEV, plasic_sz, plasic_gz);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(
            PLASIC_NLEV, plasic_sq, plasic_gq);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(1, plasic_sp, plasic_gp);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid_gradient(
            plasic_sp, plasic_gpj,
            runtime_fields->zonal_log_surface_pressure_gradient);
    }
    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error("adiabatic spherical synthesis failed (status %d)", status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    return PLASIC_RUNTIME_OK;
}

/*
 * 从完整的全局谱数组中提取当前 MPI rank 负责的局部谱块。
 *
 * 当 PLASIC_NPRO 不整除 PLASIC_NRSP 时，最后一个 rank 只有部分槽位有效：
 * 有效部分复制到局部数组，剩余槽位清零，作为不参与任何物理解算的 padding。
 *
 * Extract the local spectral block owned by the current MPI rank from the
 * complete global spectral array.
 *
 * When PLASIC_NPRO does not divide PLASIC_NRSP, the last rank owns only a
 * partial block: the valid part is copied into the local array and the
 * remaining slots are zeroed as padding that enters no physical computation.
 */
static void extract_local_spectral_block(
    const float *source_spectrum, float *destination_spectrum,
    int32_t level_count)
{
    /* 当前 rank 在每一层全局谱数组中的连续区间起点与有效长度。 */
    /* Start and valid length of this rank's contiguous range in the global
     * spectral array on each level. */
    const int32_t rank = mp_rank();
    const size_t rank_offset = (size_t)rank * (size_t)PLASIC_NSPP;
    const int32_t valid_count = plasic_local_spectral_count(rank);
    const size_t padding_count = (size_t)(PLASIC_NSPP - valid_count);
    int32_t level;

    for (level = 0; level < level_count; ++level)
    {
        float *destination =
            destination_spectrum + (size_t)level * (size_t)PLASIC_NSPP;

        if (valid_count > 0)
        {
            memcpy(
                destination,
                source_spectrum +
                    (size_t)level * (size_t)PLASIC_NRSP + rank_offset,
                (size_t)valid_count * sizeof(*destination_spectrum));
        }
        if (padding_count > 0)
        {
            memset(
                destination + valid_count, 0,
                padding_count * sizeof(*destination_spectrum));
        }
    }
}

int runtime_gridpoint_adiabatic(void)
{
    /*
     * 完成一次绝热动力的格点计算阶段。
     *
     * Perform one adiabatic dynamics.
     */
    int status;

    status = synthesize_adiabatic_fields();
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    /*
     * plasic_gp 存储的是对数气压比 π=ln(ps/psurf)
     *
     * plasic_gp stores the log pressure ratio π=ln(ps/psurf)
     */
    exp_array(
        PLASIC_NHOR, plasic_gp,
        runtime_fields->surface_pressure_ratio);

    /*
     * 在每个格点的一根 sigma 气柱内计算绝热非线性项
     *
     * Compute the adiabatic nonlinear terms within one sigma column at each
     * grid point
     */
    calcgp(
        PLASIC_NHOR, PLASIC_NLEV, plasic_kap, plasic_cpv_cpd_minus1,
        plasic_rdbrv, plasic_psurf, plasic_ww,
        plasic_rcsq, plasic_gu,
        runtime_fields->zonal_log_surface_pressure_gradient,
        plasic_gv, plasic_gpj, plasic_gq,
        runtime_fields->surface_pressure_ratio,
        plasic_gt,
        runtime_fields->nondimensional_reference_temperature,
        plasic_gd, plasic_dsigma,
        plasic_sigmah, plasic_gz, plasic_rdsig,
        plasic_t01s2, plasic_c, plasic_tkp, plasic_g,
        plasic_dw, runtime_fields->adiabatic_temperature_source,
        runtime_fields->humidity_vertical_source,
        runtime_fields->zonal_momentum_flux,
        runtime_fields->meridional_momentum_flux,
        runtime_fields->column_pressure_advection,
        runtime_fields->scaled_moist_geopotential);

    /*
     * 计算格点空间中的非线性通量乘积
     *
     * Compute nonlinear flux products in grid space
     */
    gridpointa_products(
        PLASIC_NHOR * PLASIC_NLEV, plasic_gu, plasic_gv,
        plasic_gt, plasic_gq,
        runtime_fields->scaled_moist_geopotential,
        runtime_fields->zonal_temperature_flux,
        runtime_fields->meridional_temperature_flux,
        runtime_fields->scaled_kinetic_plus_geopotential,
        runtime_fields->zonal_humidity_flux,
        runtime_fields->meridional_humidity_flux);

    /*
     * 先将显式地表气压平流 Abar_exp 从格点空间分析到谱空间。
     * 这是单层标量场，因此 level 参数为 1；分析结果先保存在完整全局谱数组中。
     *
     * First analyze the explicit surface pressure advection Abar_exp from grid
     * space to spectral space. This is a single-level scalar field, so the
     * level argument is 1; the analysis result is first stored in the complete
     * global spectral array.
     */
    status = sht_grid_to_scalar(
        1, runtime_fields->column_pressure_advection,
        runtime_fields->column_pressure_advection_spectrum);
    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error("surface-pressure tendency analysis failed (status %d)", status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }

    /* 时间推进只需要当前 rank 的谱块，因此从完整谱中提取局部的压力 tendency。 */
    /* Time stepping only needs the current rank's spectral block, so the local pressure tendency is extracted from the complete spectrum. */
    extract_local_spectral_block(
        runtime_fields->column_pressure_advection_spectrum,
        plasic_spt, 1);

    status = sht_mktend(
        PLASIC_NLEV,
        runtime_fields->adiabatic_temperature_source,
        runtime_fields->zonal_momentum_flux,
        runtime_fields->meridional_momentum_flux,
        runtime_fields->scaled_kinetic_plus_geopotential,
        runtime_fields->zonal_temperature_flux,
        runtime_fields->meridional_temperature_flux,
        runtime_fields->divergence_tendency_spectrum,
        runtime_fields->temperature_tendency_spectrum,
        runtime_fields->vorticity_tendency_spectrum);

    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error("momentum and temperature tendency analysis failed "
                          "(status %d)",
                          status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }

    status = sht_qtend(
        PLASIC_NLEV, runtime_fields->humidity_vertical_source,
        runtime_fields->zonal_humidity_flux,
        runtime_fields->meridional_humidity_flux,
        runtime_fields->humidity_tendency_spectrum);
    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error(
            "humidity tendency analysis failed (status %d)",
            status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }

    extract_local_spectral_block(
        runtime_fields->divergence_tendency_spectrum,
        plasic_sdt, PLASIC_NLEV);
    extract_local_spectral_block(
        runtime_fields->temperature_tendency_spectrum,
        plasic_stt, PLASIC_NLEV);
    extract_local_spectral_block(
        runtime_fields->vorticity_tendency_spectrum,
        plasic_szt, PLASIC_NLEV);
    extract_local_spectral_block(
        runtime_fields->humidity_tendency_spectrum,
        plasic_sqt, PLASIC_NLEV);

    gridpointa_snapshots(
        PLASIC_NHOR, PLASIC_NLEV, plasic_psurf, plasic_cv,
        plasic_rcsq, runtime_fields->surface_pressure_ratio,
        plasic_gu,
        plasic_gv,
        plasic_dp0, plasic_du0, plasic_dv0);
    return PLASIC_RUNTIME_OK;
}

static int gather_current_spectral_state(void)
{
    /*
     * 汇总时间推进后的“当前谱状态”。
     *
     * Gather the "current spectral state" after time stepping.
     */

    int status;

    /* 汇总所有垂直层上的散度谱系数，使每个 rank 都得到完整的 plasic_sd。 */
    /* Gather divergence spectral coefficients over all vertical levels so every rank obtains the complete plasic_sd. */
    status = (int)mp_allgather_spectral(
        plasic_sd, plasic_sdp, PLASIC_NRSP, PLASIC_NSPP,
        PLASIC_NLEV);
    // int32_t mp_allgather_spectral(
    //     float *global, const float *local, int32_t global_rows,
    //     int32_t local_rows, int32_t columns)

    if (status == 0)
    {
        status = (int)mp_allgather_spectral(
            plasic_sz, plasic_szp, PLASIC_NRSP, PLASIC_NSPP,
            PLASIC_NLEV);
    }
    if (status == 0)
    {
        status = (int)mp_allgather_spectral(
            plasic_st, plasic_stp, PLASIC_NRSP, PLASIC_NSPP,
            PLASIC_NLEV);
    }
    if (status == 0)
    {
        status = (int)mp_allgather_spectral(
            plasic_sp, plasic_spp, PLASIC_NRSP, PLASIC_NSPP, 1);
    }
    if (status == 0)
    {
        status = (int)mp_allgather_spectral(
            plasic_sq, plasic_sqp, PLASIC_NRSP, PLASIC_NSPP,
            PLASIC_NLEV);
    }
    if (status != 0)
    {
        runtime_set_error(
            "MPI spectral all-gather failed (status %d)",
            status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    return PLASIC_RUNTIME_OK;
}

int runtime_spectral_adiabatic(void)
{
    /*
     * 在谱空间完成一次“绝热动力”时间推进。
     *
     * Perform one "adiabatic dynamics" time step in spectral space.
     */

    /*
     * 阶段 1：求半隐式耦合方程。
     *
     * 对每个总波数 n，核心是求一根垂直气柱上的线性方程：
     *
     *   M_n Dbar = b_n，
     *
     * 其中 Dbar 是时间中心散度，plasic_bm1 保存预先计算好的 M_n^{-1}。
     * 求出 Dbar 后，再补齐整柱质量散度和完整温度 tendency。
     *
     * Stage 1: solve the semi-implicit coupled equations.
     *
     * For each total wavenumber n, the core task is to solve a linear equation
     * on one vertical column:
     *
     *   M_n Dbar = b_n,
     *
     * where Dbar is the time-centered divergence and plasic_bm1 stores the
     * precomputed M_n^{-1}. Once Dbar is obtained, the column mass divergence
     * and the complete temperature tendency are filled in.
     */
    spectrala_implicit(
        PLASIC_NSPP, PLASIC_NLEV, plasic_nadv, mp_rank(), 0,
        plasic_delt,
        plasic_nindex + mp_rank() * PLASIC_NSPP,
        plasic_g,
        runtime_fields->nondimensional_reference_temperature,
        plasic_bm1,
        plasic_sop,
        plasic_dsigma, plasic_tau, plasic_spt,
        plasic_sdt, plasic_stt,
        runtime_fields->column_mass_divergence,
        runtime_fields->centered_divergence,
        runtime_fields->complete_temperature_tendency,
        plasic_spm, plasic_sdm, plasic_szm,
        plasic_stm, plasic_sqm,
        runtime_fields->saved_surface_log_pressure,
        runtime_fields->saved_divergence,
        runtime_fields->saved_vorticity,
        runtime_fields->saved_temperature,
        runtime_fields->saved_humidity);

    /*
     * 阶段 2：完成跃蛙更新。
     *
     * 普通变量采用
     *
     *   X^(j+1) = X^(j-1) + 2Δτ R_X^j，
     *
     * 散度则由半隐式得到的时间中心值更新：
     *
     *   D^(j+1) = 2 Dbar - D^(j-1)。
     *
     * Stage 2: complete the leapfrog update.
     *
     * Ordinary variables use
     *
     *   X^(j+1) = X^(j-1) + 2Δτ R_X^j,
     *
     * while divergence is updated from the semi-implicit time-centered value:
     *
     *   D^(j+1) = 2 Dbar - D^(j-1).
     */
    spectrala_finalize(
        PLASIC_NSPP, PLASIC_NLEV, plasic_nadv, mp_rank(), 0,
        plasic_nkits, plasic_delt2, plasic_pnu,
        plasic_pnu21, runtime_fields->column_mass_divergence,
        runtime_fields->centered_divergence, plasic_szt,
        runtime_fields->complete_temperature_tendency, plasic_sqt,
        runtime_fields->saved_surface_log_pressure,
        runtime_fields->saved_divergence,
        runtime_fields->saved_vorticity,
        runtime_fields->saved_temperature,
        runtime_fields->saved_humidity,
        plasic_spm, plasic_sdm, plasic_szm,
        plasic_stm, plasic_sqm, plasic_spp,
        plasic_sdp, plasic_szp, plasic_stp,
        plasic_sqp);

    return gather_current_spectral_state();
}

static int synthesize_diabatic_fields(void)
{
    int status = sht_vortdiv_to_wind(
        PLASIC_NLEV, plasic_plavor,
        plasic_sd, plasic_sz, plasic_gu, plasic_gv);
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(
            PLASIC_NLEV, plasic_sq, plasic_gq);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(
            PLASIC_NLEV, plasic_st, plasic_gt);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_scalar_to_grid(1, plasic_sp, plasic_gp);
    }
    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error("diabatic spherical synthesis failed (status %d)", status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    return PLASIC_RUNTIME_OK;
}

int runtime_gridpoint_diabatic(int dynamics_only)
{
    /*
     * 在格点空间计算非绝热物理过程，并把结果转换回谱空间。
     *
     * Compute the diabatic physics processes in grid space and transform the
     * results back to spectral space.
     */
    int status;

    /*
     * 阶段 1：清空本步的 tendency 缓冲区，防止残留上一个时间步的物理增量。
     *
     * Stage 1: clear this step's tendency buffers to prevent residual physics
     * increments from the previous time step.
     */
    gridpointd_zero(
        PLASIC_NHOR, PLASIC_NLEV, plasic_gudt, plasic_gvdt,
        plasic_gtdt, plasic_gqdt, plasic_dudt,
        plasic_dvdt, plasic_dtdt, plasic_dqdt);

    /*
     * 阶段 2：把绝热推进后的完整谱状态合成为格点场：
     *
     * Stage 2: synthesize the complete spectral state after the adiabatic
     * advance into grid fields:
     */
    status = synthesize_diabatic_fields();
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    /*
     * 阶段 3：恢复物理参数化需要的有量纲状态：
     *
     * Stage 3: recover the dimensional state needed by the physics
     * parameterizations:
     */
    gridpointd_physical(
        PLASIC_NHOR, PLASIC_NLEV, plasic_psurf, plasic_cv,
        plasic_ct, plasic_rcsq, plasic_gp, plasic_gu,
        plasic_gv, plasic_gt,
        runtime_fields->nondimensional_reference_temperature,
        plasic_gq,
        plasic_dp, plasic_du, plasic_dv, plasic_dt,
        plasic_dq);

    /*
     * 阶段 4：依次执行启用的物理参数化。各物理过程读取上面的有量纲状态，
     * 并把产生的每秒变化率累加到 plasic_dudt/dvdt/dtdt/dqdt。
     * dynamics_only 模式下保持这些数组为零。
     *
     * Stage 4: execute the enabled physics parameterizations in sequence.
     * Each physics process reads the dimensional state above and accumulates
     * its per-second rate of change into plasic_dudt/dvdt/dtdt/dqdt. In
     * dynamics_only mode these arrays remain zero.
     */
    if (!dynamics_only)
    {
        status = plasic_physics_step();
        if (status != PLASIC_RUNTIME_OK)
        {
            runtime_set_error("an atmospheric physics kernel failed");
            return status;
        }
    }

    /*
     * 阶段 5：将有量纲物理 tendency 换回动力核的无量纲变量。
     *
     * Stage 5: convert the dimensional physics tendencies back to the
     * dynamical core's non-dimensional variables.
     */
    gridpointd_tendencies(
        PLASIC_NHOR, PLASIC_NLEV, plasic_psurf, plasic_cv,
        plasic_ct, plasic_ww, plasic_rcsq, plasic_dp,
        plasic_dudt, plasic_dvdt, plasic_dtdt,
        plasic_dqdt, plasic_gudt, plasic_gvdt,
        plasic_gtdt, plasic_gqdt);

    /*
     * 阶段 6：对格点 tendency 做球谐分析，得到完整的全球谱系数。
     *
     * Stage 6: perform spherical harmonic analysis on the grid tendencies to
     * obtain the complete global spectral coefficients.
     */
    status = sht_grid_to_scalar(
        PLASIC_NLEV, plasic_gtdt,
        runtime_fields->temperature_tendency_spectrum);
    if (status == PLASIC_SHT_OK)
    {
        status = sht_grid_to_scalar(
            PLASIC_NLEV, plasic_gqdt,
            runtime_fields->humidity_tendency_spectrum);
    }
    if (status == PLASIC_SHT_OK)
    {
        status = sht_wind_to_vortdiv(
            PLASIC_NLEV, plasic_gudt, plasic_gvdt,
            runtime_fields->divergence_tendency_spectrum,
            runtime_fields->vorticity_tendency_spectrum);
    }
    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error(
            "diabatic spherical analysis failed (status %d)",
            status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }

    extract_local_spectral_block(
        runtime_fields->temperature_tendency_spectrum,
        plasic_stt, PLASIC_NLEV);
    extract_local_spectral_block(
        runtime_fields->divergence_tendency_spectrum,
        plasic_sdt, PLASIC_NLEV);
    extract_local_spectral_block(
        runtime_fields->vorticity_tendency_spectrum,
        plasic_szt, PLASIC_NLEV);
    extract_local_spectral_block(
        runtime_fields->humidity_tendency_spectrum,
        plasic_sqt, PLASIC_NLEV);
    return PLASIC_RUNTIME_OK;
}

int runtime_spectral_diabatic(void)
{
    /*
     * 在谱空间完成非绝热物理和数值耗散，得到本时间步的最终大气状态。
     *
     * Complete the diabatic physics and numerical dissipation in spectral space
     * to obtain the final atmospheric state of this time step.
     */

    /*
     * 阶段 1：显式加入格点物理过程产生的谱 tendency。
     *
     * Stage 1: explicitly add the spectral tendencies produced by grid-point
     * physics.
     */
    spectrald_add_physics(
        PLASIC_NSPP * PLASIC_NLEV, plasic_delt2,
        plasic_szt, plasic_stt, plasic_sdt,
        plasic_sqt, plasic_szp, plasic_stp,
        plasic_sdp, plasic_sqp);

    /*
     * 阶段 2：根据加入物理过程后的候选状态，计算谱阻尼 tendency。
     *
     * Stage 2: compute the spectral damping tendencies from the candidate
     * state after adding physics.
     */
    spectrald_dissipation(
        PLASIC_NSPP, PLASIC_NLEV, mp_rank(), 0,
        plasic_nsponge,
        plasic_delt2, plasic_plavor, plasic_dampsp,
        plasic_damp, runtime_fields->temperature_diffusion_rate,
        runtime_fields->humidity_diffusion_rate,
        runtime_fields->friction_rate,
        runtime_fields->divergence_diffusion_rate,
        runtime_fields->vorticity_diffusion_rate,
        plasic_Lnkpp, plasic_srp, plasic_sdp,
        plasic_szp, plasic_stp, plasic_sqp,
        plasic_spp, plasic_stt, plasic_sqt,
        plasic_sdt, plasic_szt,
        runtime_fields->divergence_friction_tendency,
        runtime_fields->divergence_diffusion_tendency,
        runtime_fields->vorticity_friction_tendency,
        runtime_fields->vorticity_diffusion_tendency);

    spectrald_finalize(
        PLASIC_NSPP, PLASIC_NLEV, plasic_nkits,
        plasic_delt2, plasic_pnu, plasic_sdt,
        plasic_szt, plasic_stt, plasic_sqt,
        plasic_spm,
        plasic_sdm, plasic_szm, plasic_stm,
        plasic_sqm, plasic_spp, plasic_sdp,
        plasic_szp, plasic_stp, plasic_sqp);

    return gather_current_spectral_state();
}
