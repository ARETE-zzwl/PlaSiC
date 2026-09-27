// ============================================================================
// 热力海冰数值核。
// ============================================================================
#include "icemod_kernels.h"

#include <math.h>
#include <stddef.h>

// 月气候态的布局为 field[month][horizontal]；同一个月的 nhor 个格点连续。
// 用 size_t 完成乘法，避免大网格下 int32_t 索引溢出。
static size_t monthly_index(int32_t horizontal, int32_t month, int32_t nhor)
{
    return (size_t)horizontal + (size_t)month * (size_t)nhor;
}

// ----------------------------------------------------------------------------
// 汇总本步可改变海冰相态的能量，并按潜热守恒更新 h_i。
//
// 进入本函数时：
//   conductive_flux      已是 Q_c^raw，加上融雪后的剩余能量；
//   freezing_flux        是 Q_f^adj = K_z (T_f - T_o)；
//   ocean_flux           是海洋冻结点约束返回的 Q_o->i；
//   snow_conversion_flux 是雪转冰对应的 Q_s->i。
//
// 因而此处构造 Q_ice = conductive_flux + snow_conversion_flux
//                           - freezing_flux + ocean_flux。
// 正 Q_ice 使冰变薄，负 Q_ice 使冰增长。
// 若 h_i 试算为负，超额融化，能量会显式回补到 conductive_flux / negative_ice_flux，而不是被截断丢失。
// ----------------------------------------------------------------------------
void ice_make_ice(
    int32_t nhor, float tfreeze, float ice_density, float latent_heat_ice,
    float timestep, const float *land_mask, const float *ocean_flux,
    const float *sst, const float *snow_conversion_flux,
    const float *freezing_flux, float *conductive_flux,
    float *ice_melt_flux, float *negative_ice_flux, float *ice_thickness)
{
    /*
     * 每 1 m 海冰在一个时间步内对应的潜热通量：
     *     ice_energy_per_timestep = rho_i L_i / dt.
     * 因此，后面用“相变通量 / ice_energy_per_timestep”就可以得到
     * 本步的冰厚变化（单位：m）。
     */
    const float ice_energy_per_timestep = ice_density * latent_heat_ice / timestep;

    /* 每个水平格点都是一个彼此独立的冰雪柱。 */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /*
         * ocean_flux 是海洋侧冻结点约束返回给冰的通量。它只存在于海洋
         * 格点，并先并入 conductive_flux，之后和其它冰相变能量统一结算。
         */
        if (land_mask[horizontal] < 0.5f)
        {
            conductive_flux[horizontal] =
                conductive_flux[horizontal] + ocean_flux[horizontal];
        }

        /*
         * 只有以下两种情况允许改变冰的潜热储量：
         *
         *   1. 海水已经达到冻结点，可以生成新冰；
         *   2. 已经存在海冰，可以融冰或继续增厚。
         *
         * 无冰且 SST 高于冻结点时，热量只能留给海洋显热，不能用于
         * 凭空生成海冰。
         */
        if (sst[horizontal] <= tfreeze || ice_thickness[horizontal] > 0.0f)
        {
            /*
             * 海冰相变总账：
             *   Q_ice = Q_c + Q_s->i - Q_f,
             */
            ice_melt_flux[horizontal] = conductive_flux[horizontal] +
                                        snow_conversion_flux[horizontal] -
                                        freezing_flux[horizontal];

            /* 潜热守恒：h_i <- h_i - Q_ice dt / (rho_i L_i)。 */
            ice_thickness[horizontal] = ice_thickness[horizontal] -
                                        ice_melt_flux[horizontal] /
                                            ice_energy_per_timestep;

            /*
             * Q_c - Q_f 用于改变冰厚，剩下的 Q_f 仍是海洋侧的接口热通量。
             * 调用方随后会把 conductive_flux 交回 ocean 模块。
             */
            conductive_flux[horizontal] = freezing_flux[horizontal];
        }

        /*
         * 如果试算冰厚不为正，说明融化能量超过了现有冰的全部潜热。
         * h_i 不能变成负数；超出的能量必须退回海洋预算，而不能丢失。
         */
        if (ice_thickness[horizontal] <= 0.0f)
        {
            /*
             * 下面三行保持原有的浮点运算和符号关系：
             *   - conductive_flux：增加退回海洋的超额能量；
             *   - ice_melt_flux：从请求融化量中扣除超额部分；
             *   - negative_ice_flux：记录超额融化能量。
             */
            conductive_flux[horizontal] = conductive_flux[horizontal] -
                                          ice_thickness[horizontal] *
                                              ice_energy_per_timestep;
            ice_melt_flux[horizontal] = ice_melt_flux[horizontal] +
                                        ice_thickness[horizontal] *
                                            ice_energy_per_timestep;
            negative_ice_flux[horizontal] = negative_ice_flux[horizontal] -
                                            ice_thickness[horizontal] *
                                                ice_energy_per_timestep;
            ice_thickness[horizontal] = 0.0f;
        }

        /* 避免浮点舍入留下会误触发“已有冰”分支的幽灵冰。 */
        if (fabsf(ice_thickness[horizontal]) < 1.0e-8f)
        {
            ice_thickness[horizontal] = 0.0f;
        }
    }
}

// ----------------------------------------------------------------------------
// 将厚度变化映射为连续紧密度 c_i^*。
// 它用于低成本地表示冰缘扩张/退缩；最终传给 atmosphere/sea 的 ice_cover
// 会在调用方按阈值二值化，不能解释为格点内的真实部分覆盖面积。
// ----------------------------------------------------------------------------
void ice_compactness(
    int32_t nhor, const float *old_thickness, const float *new_thickness,
    float *compactness)
{
    // h_r = 0.5 m：经验厚度尺度，控制新冰向未覆盖面积扩张的速度。
    const float reference_thickness = 0.5f;

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        // h_i 增加时，只让尚未覆盖的 1-c_i^* 部分扩张。
        if (new_thickness[horizontal] > old_thickness[horizontal])
        {
            compactness[horizontal] = compactness[horizontal] +
                                      (1.0f - compactness[horizontal]) *
                                          (new_thickness[horizontal] -
                                           old_thickness[horizontal]) /
                                          reference_thickness;
            compactness[horizontal] =
                fminf(compactness[horizontal], 1.0f);
        }

        // h_i 减少时，c_i^* 按相对厚度损失衰减；此分支只有 old_thickness > 0
        // 时才有物理意义，调用顺序保证该条件。
        if (new_thickness[horizontal] < old_thickness[horizontal])
        {
            compactness[horizontal] = compactness[horizontal] +
                                      compactness[horizontal] *
                                          (new_thickness[horizontal] -
                                           old_thickness[horizontal]) /
                                          (2.0f * old_thickness[horizontal]);
            compactness[horizontal] =
                fmaxf(compactness[horizontal], 0.0f);
        }
        // 冰的潜热储库耗尽时，连续紧密度也必须归零。
        if (new_thickness[horizontal] <= 0.0f)
        {
            compactness[horizontal] = 0.0f;
        }
    }
}

// ----------------------------------------------------------------------------
// 更新冰上雪水当量 S，并执行“先融雪、再融冰”和雪转冰两条能量路径。
//
// ice_melt_energy 是冰表达到 T_m 后不能继续升温的剩余能量 Q_m。
// 本函数首先以它融雪；只有积雪耗尽后的余量才加入 conductive_flux，供
// ice_make_ice 融冰。snow_conversion_flux 则是重雪压沉后雪转冰的 Q_s->i，
// 同样由 ice_make_ice 统一写入冰的潜热预算。
// ----------------------------------------------------------------------------
void ice_snow(
    int32_t nhor, float minimum_ice_cover, float seawater_density,
    float ice_density, float latent_heat_snow, float timestep,
    const float *land_mask, const float *ice_cover, float *snowfall,
    const float *ice_melt_energy, float *ice_thickness, float *snow,
    float *conductive_flux, float *snow_melt, float *snow_melt_flux,
    float *snow_conversion_flux, float *snow_to_ice)
{
    // 所有过程均为本地一维冰雪柱过程，格点之间没有雪或海冰输送。
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        // c_i 未达到保雪阈值时，把降雪视为落在开阔水面：不建立 S，
        // 而是立即从海面预算扣除其融化潜热。
        if (ice_cover[horizontal] < minimum_ice_cover &&
            snowfall[horizontal] > 0.0f)
        {
            // 1000 kg m^-3 将雪水当量转换为面质量；负号表示这份能量不能再用于冰的相变。
            conductive_flux[horizontal] = conductive_flux[horizontal] -
                                          snowfall[horizontal] * 1000.0f *
                                              latent_heat_snow;
            snow_melt[horizontal] = snowfall[horizontal] + snow_melt[horizontal];
            snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                         snowfall[horizontal] * 1000.0f *
                                             latent_heat_snow;
            snowfall[horizontal] = 0.0f;
        }

        // 存在冰时：Q_m 先改变 S，剩余部分才能改变 h_i。
        if (ice_thickness[horizontal] > 0.0f)
        {
            // rho_w L_s S/dt：融尽全部积雪所需的通量。
            const float snow_energy = snow[horizontal] * 1000.0f *
                                      latent_heat_snow / timestep;
            const int32_t melts_all_snow =
                ice_melt_energy[horizontal] > snow_energy;
            const int32_t melts_some_snow =
                !melts_all_snow && ice_melt_energy[horizontal] > 0.0f &&
                snow[horizontal] > 0.0f;

            // 若 Q_m 足够融尽 S，则记录未被雪消耗的余量；否则余量为零。
            const float remaining_melt_energy =
                melts_all_snow
                    ? ice_melt_energy[horizontal] -
                          snow[horizontal] * 1000.0f * latent_heat_snow / timestep
                    : (melts_some_snow ? 0.0f : ice_melt_energy[horizontal]);

            if (melts_all_snow)
            {
                snow_melt[horizontal] = snow_melt[horizontal] +
                                        snow[horizontal] / timestep;
                snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                             snow[horizontal] * 1000.0f *
                                                 latent_heat_snow / timestep;
                snow[horizontal] = 0.0f;
            }
            else if (melts_some_snow)
            {
                snow_melt[horizontal] = snow_melt[horizontal] +
                                        ice_melt_energy[horizontal] /
                                            latent_heat_snow / 1000.0f;
                snow[horizontal] = snow[horizontal] -
                                   ice_melt_energy[horizontal] * timestep /
                                       latent_heat_snow / 1000.0f;
                snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                             ice_melt_energy[horizontal];
            }
            // 仅雪融尽后，Q_m^rem 才会成为可融冰的通量。
            conductive_flux[horizontal] = conductive_flux[horizontal] +
                                          remaining_melt_energy;
        }
        // 无冰时，冰上雪水当量没有可依附的储库；将残余 S 全部移回能量预算。
        else
        {
            snow_melt[horizontal] = snow_melt[horizontal] +
                                    snow[horizontal] / timestep;
            conductive_flux[horizontal] = conductive_flux[horizontal] -
                                          snow[horizontal] * 1000.0f *
                                              latent_heat_snow / timestep;
            snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                         snow[horizontal] * 1000.0f *
                                             latent_heat_snow / timestep;
            snow[horizontal] = 0.0f;
        }

        // snow_to_ice 是质量转换诊断；每步重新计算，不能沿用旧值。
        snow_to_ice[horizontal] = 0.0f;
        // 由阿基米德浮力条件 rho_w S > (rho_o-rho_i) h_i
        // 判断积雪是否压低冰雪界面并触发雪转冰。
        if (land_mask[horizontal] < 0.5f &&
            1000.0f * snow[horizontal] >
                (seawater_density - ice_density) * ice_thickness[horizontal])
        {
            // 目标厚度 h_i,target = (rho_w S + rho_i h_i)/rho_o，
            // 即转换后冰雪界面恰回到海平面。
            const float new_ice =
                (1000.0f * snow[horizontal] +
                 ice_density * ice_thickness[horizontal]) /
                seawater_density;
            // 将相同质量从 S 移至 h_i：snow_change 通常为负，ice_change 为正。
            const float ice_change = new_ice - ice_thickness[horizontal];
            const float snow_change = -ice_change * ice_density / 1000.0f;

            // 此处只更新雪；冰的潜热变化必须和其他 Q_ice 项一起在 make_ice 中结算。
            snow[horizontal] = snow[horizontal] + snow_change;
            
            // 将雪的质量变化转为 Q_s->i。
            snow_conversion_flux[horizontal] =
                snow_change * 1000.0f * latent_heat_snow / timestep;
            snow_to_ice[horizontal] =
                ice_change * ice_density / 1000.0f / timestep;
            // 这部分雪是转冰而非融化，故从融雪诊断中扣回。
            snow_melt[horizontal] = snow_melt[horizontal] -
                                    snow_change / timestep;
            snow_melt_flux[horizontal] = snow_melt_flux[horizontal] -
                                         snow_change * 1000.0f *
                                             latent_heat_snow / timestep;
        }
    }
}

// ----------------------------------------------------------------------------
// 诊断冰雪柱的原始导热 Q_c^raw 和冻结点调整 Q_f^adj。
//
// 雪水当量先转换为几何雪深 h_s = rho_w S/rho_s，再使用
// 厚度加权导热率和 K_z = k_bar/(h_i+h_s)。于是
//
//   diagnosed_conductive_flux = Q_c^raw = K_z (T_i - T_o)
//   freezing_flux             = Q_f^adj = K_z (T_f - T_o)
//
// 在预报冰的默认路径中，ice_make_ice 取两者之差，得到由 T_i-T_f 决定的
// 导热相变贡献。薄冰不解析该冰雪柱；达到 maximum_thickness
// 时则按既有参数化关闭垂直导热。
// ----------------------------------------------------------------------------
void ice_conductive_flux(
    int32_t nhor, int32_t ocean_flux_correction, int32_t prognostic_ice,
    float minimum_thickness, float maximum_thickness, float tfreeze,
    float snow_density, float ice_conductivity, float snow_conductivity,
    const float *ice_thickness, const float *snow,
    const float *surface_temperature, const float *sst,
    const float *atmospheric_heat, float *diagnosed_conductive_flux,
    float *conductive_flux, float *freezing_flux)
{
    // 逐格点计算等效导热系数及相应通量。
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        diagnosed_conductive_flux[horizontal] = 0.0f;
        // h_i < h_min 时没有独立、可解析的冰雪柱；大气热量直接进入后续预算。
        if (ice_thickness[horizontal] < minimum_thickness)
        {
            conductive_flux[horizontal] = atmospheric_heat[horizontal];
            freezing_flux[horizontal] = 0.0f;
        }
        // 厚冰：以串联冰雪层的等效导热参数化连接 T_i 和 T_o。
        else
        {
            const float snow_depth = 1000.0f / snow_density * snow[horizontal];
            const float layer_mean_conductivity =
                (snow_conductivity * snow_depth +
                 ice_conductivity * ice_thickness[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);
            const float effective_conductivity =
                maximum_thickness >= 0.0f &&
                        ice_thickness[horizontal] >= maximum_thickness
                    ? 0.0f
                    : layer_mean_conductivity;
            diagnosed_conductive_flux[horizontal] =
                effective_conductivity *
                (surface_temperature[horizontal] - sst[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);
            conductive_flux[horizontal] = diagnosed_conductive_flux[horizontal];
            freezing_flux[horizontal] =
                effective_conductivity * (tfreeze - sst[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);
        }
    }

    // 气候态路径或海洋 SST 在线订正路径不使用冻结点调整，
    // 故显式关掉 Q_f^adj。该开关与 OCEANMOD 共享。
    if (ocean_flux_correction == 1 || prognostic_ice == 0)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            freezing_flux[horizontal] = 0.0f;
        }
    }
}


// ----------------------------------------------------------------------------
// 更新冰雪表面温度 T_i，并把到达融点后不能继续升温的能量转为 Q_m。
//
// 对 h_i >= h_min 的厚冰，表层有效热容量取
// C_i = rho_i c_p,i h_min，并半隐式求解
//
//   T_i^(n+1) = [(C_i/dt)T_i^n + Q_a + K_z T_o]
//                / [C_i/dt + K_z]               
//
// 若试算温度超过 T_m，则钳制到 T_m，并输出融雪/融冰候选通量 melt_energy。
// 薄冰不建立独立表层储热，直接令 T_i=T_o。
// ----------------------------------------------------------------------------
void ice_skin_temperature(
    int32_t nhor, float minimum_thickness, float maximum_thickness,
    float melting_temperature, float ice_density, float ice_heat_capacity,
    float snow_density, float ice_conductivity, float snow_conductivity,
    float timestep, const float *ice_thickness, const float *snow,
    const float *sst, const float *atmospheric_heat,
    float *surface_temperature, float *surface_storage_flux,
    float *melt_energy)
{
    // 零层模型只给最薄的有效表皮配置热容量，而非求解多层冰内温度剖面。
    const float layer_heat_capacity =
        minimum_thickness * ice_density * ice_heat_capacity;

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        surface_storage_flux[horizontal] = 0.0f;
        melt_energy[horizontal] = 0.0f;
        
        // T_i^n 不允许高于 T_m；超过融点的旧值没有可保留的冰雪显热意义。
        const float old_temperature =
            fminf(melting_temperature, surface_temperature[horizontal]);
        const float snow_depth = 1000.0f / snow_density * snow[horizontal];

        // 厚冰
        if (ice_thickness[horizontal] >= minimum_thickness)
        {
            const float heat_capacity_per_timestep =
                layer_heat_capacity / timestep;
            const float mean_conductivity =
                (snow_conductivity * snow_depth +
                 ice_conductivity * ice_thickness[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);

            const float unconstrained_conductivity_per_depth =
                mean_conductivity /
                (ice_thickness[horizontal] + snow_depth);

            const float effective_conductivity_per_depth =
                ice_thickness[horizontal] >= maximum_thickness &&
                        maximum_thickness >= 0.0f
                    ? 0.0f
                    : unconstrained_conductivity_per_depth;

            // 分子中的 Q_a + K_z T_o。
            const float total_flux = atmospheric_heat[horizontal] +
                                     effective_conductivity_per_depth *
                                         sst[horizontal];

            surface_temperature[horizontal] =
                (heat_capacity_per_timestep * old_temperature + total_flux) /
                (heat_capacity_per_timestep + effective_conductivity_per_depth);

            // T_i 不能高于 T_m：剩余能量按先融雪、后融冰的次序处理。
            if (surface_temperature[horizontal] > melting_temperature)
            {
                // Q_m = Q_a + K_z(T_o-T_m) - (C_i/dt)(T_m-T_i^n)，即式 (26)。
                melt_energy[horizontal] =
                    atmospheric_heat[horizontal] +
                    effective_conductivity_per_depth *
                        (sst[horizontal] - melting_temperature) -
                    heat_capacity_per_timestep *
                        (melting_temperature - old_temperature);
                surface_temperature[horizontal] = melting_temperature;
            }
            // C_i (T_i^(n+1)-T_i^n)/dt：本步留在冰雪表皮的显热通量。
            surface_storage_flux[horizontal] =
                heat_capacity_per_timestep *
                (surface_temperature[horizontal] - old_temperature);
        }
        // 薄冰/无冰：海水就是大气看到的表面，不设独立冰表温度。
        else
        {
            surface_temperature[horizontal] = sst[horizontal];
            // 虽未达到可解析厚冰阈值，残余薄冰仍可由大气正通量融尽。
            if (ice_thickness[horizontal] > 0.0f)
            {
                melt_energy[horizontal] = atmospheric_heat[horizontal];
            }
        }
    }
}

// ----------------------------------------------------------------------------
// 在两个相邻月气候态间插值。
// ----------------------------------------------------------------------------
void ice_interpolate_climatology(
    int32_t nhor, int32_t previous_month1, int32_t previous_month2,
    float previous_weight1, float previous_weight2, int32_t month1,
    int32_t month2, float weight1, float weight2,
    const float *climatological_sst, const float *climatological_cover,
    const float *climatological_thickness, float *previous_sst,
    float *current_sst, float *current_cover, float *current_thickness)
{
    /*
     * 这个函数只负责“查表并插值”，不推进海冰状态：
     *
     *   - climatological_* 是按月份保存的固定气候态数据；
     *   - month1/month2 是当前时刻两侧的月份；
     *   - weight1/weight2 是两个月份的插值权重，通常满足
     *     weight1 + weight2 = 1；
     *   - previous_* 表示前一个模型时刻使用的月份和权重。
     *
     * 月份数组的内存布局是 field[month][horizontal]
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /*
         * 保存前一个模型时刻的气候态海表温度。
         * 它主要用于需要比较相邻模型时刻的季节性过程。
         */
        previous_sst[horizontal] =
            previous_weight1 *
                climatological_sst[monthly_index(horizontal, previous_month1, nhor)] +
            previous_weight2 *
                climatological_sst[monthly_index(horizontal, previous_month2, nhor)];

        /* 当前模型时刻的气候态海表温度。 */
        current_sst[horizontal] =
            weight1 * climatological_sst[monthly_index(horizontal, month1, nhor)] +
            weight2 * climatological_sst[monthly_index(horizontal, month2, nhor)];

        /* 当前模型时刻的气候态海冰厚度。 */
        current_thickness[horizontal] =
            weight1 *
                climatological_thickness[monthly_index(horizontal, month1, nhor)] +
            weight2 *
                climatological_thickness[monthly_index(horizontal, month2, nhor)];

        /* 当前模型时刻的气候态海冰覆盖率。 */
        current_cover[horizontal] =
            weight1 * climatological_cover[monthly_index(horizontal, month1, nhor)] +
            weight2 * climatological_cover[monthly_index(horizontal, month2, nhor)];

        // 数值防护：厚度和覆盖率是非负状态量。
        current_cover[horizontal] =
            fmaxf(current_cover[horizontal], 0.0f);
        current_thickness[horizontal] =
            fmaxf(current_thickness[horizontal], 0.0f);

        /* 没有冰厚时，覆盖率也必须为零，保持两个状态量的一致性。 */
        // 没有潜热储库就不能声明存在海冰覆盖。
        if (current_thickness[horizontal] <= 0.0f)
        {
            current_cover[horizontal] = 0.0f;
        }
    }
}


// ----------------------------------------------------------------------------
// 从月气候态覆盖率构造经验海冰厚度，供非预报冰或冷启动初始化使用。
//
// 这里的 thickness 不是通过海冰潜热方程预报出来的厚度，而是根据给定的
// 连续覆盖率 compactness 诊断出的“参考厚度”。北、南两个纬向半球使用不同
// 的经验映射；北侧还乘以一个逐月修正因子，以表示季节性厚度差异。
//
// 三个二维/三维量的布局都是：
//
//     field[month][latitude][longitude]
//
// northern_factor 有 14 个元素，包含月气候态使用的年界 halo；调用方必须
// 保证 month 索引位于该数组范围内。
// ----------------------------------------------------------------------------
void ice_make_thickness(
    int32_t nlon, int32_t nlat, int32_t nmonths,
    const float *compactness, float *thickness)
{
    /* 北侧海冰厚度的逐月经验修正；最后两个值对应年界 halo 月份。 */
    static const float northern_factor[14] = {
        0.912f, 0.942f, 1.0f, 1.058f, 1.124f, 1.161f, 1.175f,
        1.058f, 0.931f, 0.883f, 0.88f, 0.876f, 0.912f, 0.942f};
    /*
     * 北侧经验关系：
     *   cover < 0.1       ：只有少量冰，赋予最小厚度；
     *   0.1 <= cover < 0.9：在最小和最大厚度之间线性插值；
     *   cover >= 0.9      ：赋予最大厚度。
     */
    const float northern_minimum_cover = 0.1f;
    const float northern_maximum_cover = 0.9f;
    const float northern_minimum_thickness = 0.25f;
    const float northern_maximum_thickness = 3.0f;
    /*
     * 南侧经验关系：
     *   cover < 0.25       ：厚度为 0.25 + cover；
     *   cover >= 0.25      ：厚度固定为 0.50 m。
     * 南侧不使用月修正。
     */
    const float southern_minimum_cover = 0.25f;
    const float southern_minimum_thickness = 0.25f;
    const float southern_maximum_thickness = 0.50f;

    /* 一个“月份 × 纬度 × 经度”平面的格点数。 */
    const size_t plane = (size_t)nlon * (size_t)nlat;

    /*
     * 先清零全部输出。这样 cover <= 0 的格点无需显式赋值，仍保持无冰；
     * 同时也保证输出数组中不会残留调用前的旧厚度。
     */
    for (size_t index = 0; index < plane * (size_t)nmonths; ++index)
    {
        thickness[index] = 0.0f;
    }

    /* 按 month → latitude → longitude 遍历，最内层 longitude 与内存布局一致。 */
    for (int32_t month = 0; month < nmonths; ++month)
    {
        for (int32_t latitude = 0; latitude < nlat; ++latitude)
        {
            for (int32_t longitude = 0; longitude < nlon; ++longitude)
            {
                /* 将三维下标转换为 field[month][latitude][longitude] 的一维下标。 */
                const size_t index = (size_t)longitude +
                                     (size_t)latitude * (size_t)nlon +
                                     (size_t)month * plane;
                const float cover = compactness[index];
                
                /* 数组前半纬向带采用北侧经验关系；真实南北方向由网格纬度排序决定。 */
                if (latitude < nlat / 2)
                {
                    /* 覆盖率很高时直接采用北侧最大厚度，再乘季节因子。 */
                    if (cover >= northern_maximum_cover)
                    {
                        thickness[index] =
                            northern_maximum_thickness * northern_factor[month];
                    }

                    /* 中等覆盖率：在最小厚度和最大厚度之间线性插值。 */
                    else if (cover >= northern_minimum_cover)
                    {
                        thickness[index] =
                            (northern_minimum_thickness +
                             (northern_maximum_thickness -
                              northern_minimum_thickness) *
                                 (cover - northern_minimum_cover) /
                                 (northern_maximum_cover -
                                  northern_minimum_cover)) *
                            northern_factor[month];
                    }
                    
                    /* 只要存在少量海冰，就至少赋予北侧最小厚度。 */
                    else if (cover > 0.0f)
                    {
                        thickness[index] =
                            northern_minimum_thickness * northern_factor[month];
                    }
                }
                /* 后半纬向带采用较薄、无逐月因子的南侧经验关系。 */
                else
                {
                    /* 南侧覆盖率达到阈值后使用固定的最大经验厚度。 */
                    if (cover >= southern_minimum_cover)
                    {
                        thickness[index] = southern_maximum_thickness;
                    }
                    /* 南侧低覆盖率时，厚度按 0.25 + cover 增加。 */
                    else if (cover > 0.0f)
                    {
                        thickness[index] =
                            southern_minimum_thickness + cover;
                    }
                }
            }
        }
    }
}
