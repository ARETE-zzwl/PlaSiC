#include "plasic_kernels.h"
#include "math_tools/linear_algebra.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

/*
 * 本文件实现 PlaSiC 动力核中不依赖球谐库的数值算子。阅读时可按下列主线理解：
 *
 *   谱状态 -> 格点非线性动力(calcgp 位于同目录的 calgp.c)
 *           -> 显式谱 tendency
 *           -> 半隐式重力波修正(spectrala_*)
 *           -> 跃蛙更新和 Robert-Asselin 滤波
 *           -> 物理 tendency、谱耗散与最终状态(spectrald_*)。
 *
 * 水平场使用两类表示：谱数组的 point 实际是谱系数下标，格点数组的 point 是
 * 本 MPI 分块内的高斯格点下标；二者都采用 [level][point] 的层优先布局。
 * vertical_g、vertical_c、vertical_tau 是半隐式求解中使用的垂直算子。
 * 同一数组在整个生命周期内只表示一种状态或物理量；正常的时间推进与同类
 * tendency 累加仍写回该状态，但不会再借用它存放另一种中间量。
 */

static size_t field_index(int32_t point, int32_t level, int32_t point_count)
{
    return (size_t)level * (size_t)point_count + (size_t)point;
}

static float integer_power(float base, int32_t exponent);

// 从 π=ln(ps/p0) 单独生成 exp(π)=ps/p0；输入始终保留对数气压语义。
void exp_array(
    int32_t count, const float *log_values, float *exponential_values)
{
    int32_t index;
    for (index = 0; index < count; ++index)
    {
        exponential_values[index] = expf(log_values[index]);
    }
}

/*
 * 在格点空间形成非线性乘积，再交给球谐分析计算水平通量散度/旋度：
 *   U T'、V T' 用于温度平流，Uχ、Vχ 用于质量加权水汽输送；
 *   scaled_kinetic_plus_geopotential 保存
 *   U²+V²+2cos²φ·δΦ = 2cos²φ(K+δΦ)。
 * 这里故意不乘 1/2：公共的 2cos²φ 缩放由后续谱适配层统一消去。
 */
void gridpointa_products(
    int32_t count, const float *zonal_wind, const float *meridional_wind,
    const float *temperature, const float *humidity,
    const float *scaled_moist_geopotential, float *zonal_temperature_flux,
    float *meridional_temperature_flux,
    float *scaled_kinetic_plus_geopotential,
    float *zonal_humidity_flux, float *meridional_humidity_flux)
{
    int32_t index;
    for (index = 0; index < count; ++index)
    {
        zonal_temperature_flux[index] = zonal_wind[index] * temperature[index];
        meridional_temperature_flux[index] = meridional_wind[index] * temperature[index];

        const float robert_wind_squared =
            zonal_wind[index] * zonal_wind[index] +
            meridional_wind[index] * meridional_wind[index];

        zonal_humidity_flux[index] = zonal_wind[index] * humidity[index];
        meridional_humidity_flux[index] = meridional_wind[index] * humidity[index];

        scaled_kinetic_plus_geopotential[index] =
            robert_wind_squared + scaled_moist_geopotential[index];
    }
}

/*
 * 保存进入物理过程前的有量纲状态快照。surface_pressure 已是 ps/p0，
 * Robert 风乘 sqrt(rcsq)=1/cosφ 后还原物理风，再乘各自的量纲尺度。
 * 快照用于比较物理过程前后的状态或诊断耗散能量。
 */
void gridpointa_snapshots(
    int32_t horizontal_count, int32_t level_count, float reference_pressure,
    float velocity_scale, const float *rcsq,
    const float *surface_pressure_ratio,
    const float *zonal_wind, const float *meridional_wind,
    float *saved_surface_pressure, float *saved_zonal_wind,
    float *saved_meridional_wind)
{
    int32_t point;
    int32_t level;
    for (point = 0; point < horizontal_count; ++point)
    {
        saved_surface_pressure[point] = reference_pressure * surface_pressure_ratio[point];
    }

    for (level = 0; level < level_count; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index = field_index(point, level, horizontal_count);
            const float metric = sqrtf(rcsq[point]);
            saved_zonal_wind[index] =
                velocity_scale * zonal_wind[index] * metric;
            saved_meridional_wind[index] =
                velocity_scale * meridional_wind[index] * metric;
        }
    }
}

// 简单线性求和

// 从数组每个元素减去同一常数，常用于移除某个全局平均/基准模态。

/*
 * 清零格点非绝热步骤的各类 tendency 缓冲区。动力数组只有 L 个整层；物理数组
 * 按历史接口分配 L+1 个平面（额外平面可保存地表），因此两段循环长度不同。
 */
void gridpointd_zero(
    int32_t horizontal_count, int32_t level_count,
    float *zonal_tendency, float *meridional_tendency,
    float *temperature_tendency, float *humidity_tendency,
    float *physical_zonal_tendency, float *physical_meridional_tendency,
    float *physical_temperature_tendency, float *physical_humidity_tendency)
{
    const int32_t count = horizontal_count * level_count;
    const int32_t physical_count = horizontal_count * (level_count + 1);
    int32_t index;

    for (index = 0; index < count; ++index)
    {
        zonal_tendency[index] = 0.0f;
        meridional_tendency[index] = 0.0f;
        temperature_tendency[index] = 0.0f;
        humidity_tendency[index] = 0.0f;
    }

    for (index = 0; index < physical_count; ++index)
    {
        physical_zonal_tendency[index] = 0.0f;
        physical_meridional_tendency[index] = 0.0f;
        physical_temperature_tendency[index] = 0.0f;
        physical_humidity_tendency[index] = 0.0f;
    }
}

/*
 * 把动力核无量纲格点状态转换成物理参数化使用的有量纲量：
 *   ps=p0 expπ，(u,v)=velocity_scale·(U,V)/cosφ，T=temperature_scale(T'+T0)，
 *   q=χ p0/ps。
 */
void gridpointd_physical(
    int32_t horizontal_count, int32_t level_count,
    float reference_pressure, float velocity_scale, float temperature_scale,
    const float *rcsq, const float *log_surface_pressure,
    const float *zonal_wind, const float *meridional_wind,
    const float *temperature, const float *reference_temperature,
    const float *humidity, float *surface_pressure,
    float *physical_zonal_wind, float *physical_meridional_wind,
    float *physical_temperature, float *physical_humidity)
{
    int32_t point;
    int32_t level;
    for (point = 0; point < horizontal_count; ++point)
    {
        surface_pressure[point] =
            reference_pressure * expf(log_surface_pressure[point]);
    }

    for (level = 0; level < level_count; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index = field_index(point, level, horizontal_count);
            const float metric = sqrtf(rcsq[point]);

            physical_zonal_wind[index] =
                velocity_scale * zonal_wind[index] * metric;
            physical_meridional_wind[index] =
                velocity_scale * meridional_wind[index] * metric;

            physical_temperature[index] = temperature_scale * (temperature[index] + reference_temperature[level]);

            physical_humidity[index] =
                humidity[index] * reference_pressure / surface_pressure[point];
        }
    }
}

/*
 * 将物理参数化产生的有量纲 tendency 换回动力核的无量纲变量并累加：
 * 风需去掉 Robert 度量和速度/时间尺度，温度去掉温度/时间尺度；水汽则由
 * dq/dt 转为 dχ/dτ=(ps/p0)dq/dτ。这里使用 +=，保留此前形成的绝热 tendency。
 */
void gridpointd_tendencies(
    int32_t horizontal_count, int32_t level_count,
    float reference_pressure, float velocity_scale, float temperature_scale,
    float rotation_rate, const float *rcsq, const float *surface_pressure,
    const float *physical_zonal_tendency,
    const float *physical_meridional_tendency,
    const float *physical_temperature_tendency,
    const float *physical_humidity_tendency,
    float *zonal_tendency, float *meridional_tendency,
    float *temperature_tendency, float *humidity_tendency)
{
    int32_t level;
    int32_t point;
    for (level = 0; level < level_count; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index = field_index(point, level, horizontal_count);
            const float metric = sqrtf(rcsq[point]);

            zonal_tendency[index] =
                physical_zonal_tendency[index] / metric / velocity_scale / rotation_rate + zonal_tendency[index];

            meridional_tendency[index] = physical_meridional_tendency[index] / metric / velocity_scale / rotation_rate + meridional_tendency[index];

            temperature_tendency[index] = physical_temperature_tendency[index] / temperature_scale / rotation_rate + temperature_tendency[index];

            humidity_tendency[index] = physical_humidity_tendency[index] * surface_pressure[point] / rotation_rate / reference_pressure + humidity_tendency[index];
        }
    }
}

/*
 * 求半隐式重力波方程，并保存跃蛙更新所需的旧层。
 *
 * 对每个 n>0 的球谐模态，垂直中心散度 Dbar 满足
 *   M_n Dbar = b_n,
 *   M_n = I/[n(n+1)] + Δτ²(G·Tau + T0·Δσ^T)。
 * inverse_matrix 已由 makebm 按总波数 n 预先求逆；同一 n 的所有 m 共用它。
 * 本函数组装 b_n、做一次 L×L 矩阵向量乘法，并用求得的 Dbar 补齐
 *   Abar = Abar_exp + Σ_k Δσ_k Dbar_k，
 *   R_T = R_T_exp - Tau·Dbar。
 * 显式 tendency、完整柱质量散度、完整温度 tendency 和时间中心散度分别使用
 * 独立数组，避免任何数组在函数执行中改变物理含义。
 */
void spectrala_implicit(
    int32_t spectral_count, int32_t level_count, int32_t advection_enabled,
    int32_t process_id, int32_t root_process, float time_step,
    const int32_t *wavenumber, const float *vertical_g,
    const float *reference_temperature, const float *inverse_matrix,
    const float *surface_geopotential, const float *sigma_thickness,
    const float *vertical_tau, const float *column_pressure_advection,
    const float *explicit_divergence_tendency,
    const float *explicit_temperature_tendency,
    float *column_mass_divergence, float *centered_divergence,
    float *complete_temperature_tendency,
    const float *old_surface_log_pressure, const float *old_divergence,
    const float *old_vorticity, const float *old_temperature,
    const float *old_humidity, float *saved_surface_log_pressure,
    float *saved_divergence, float *saved_vorticity,
    float *saved_temperature, float *saved_humidity)
{
    const int32_t level_values = spectral_count * level_count;

    // 两个临时列都按 [level][spectral mode] 展平：
    //   g_times_tendency   = G·R_T'^exp，
    //   g_times_temperature = G·T'^{j-1}。
    // 它们都用于构造右端 b_n
    float *g_times_tendency =
        (float *)malloc((size_t)level_values * sizeof(float));
    float *g_times_temperature =
        (float *)malloc((size_t)level_values * sizeof(float));

    int32_t index;
    int32_t level;

    // 阶段 0：保存 j-1 时间层。
    // 本函数随后只计算 Dbar 与完整 tendency；真正的 X^{j+1} 更新由
    // spectrala_finalize 完成。因此此处先保存所有旧状态，防止调用者的状态数组
    // 在后续步骤被改写后丢失跃蛙的基点 X^{j-1}。
    // 显式 tendency 参数均为 const，语义上始终表示当前层 j 的显式右端。
    for (index = 0; index < spectral_count; ++index)
    {
        saved_surface_log_pressure[index] =
            old_surface_log_pressure[index];
    }

    for (index = 0; index < level_values; ++index)
    {
        saved_divergence[index] = old_divergence[index];
        saved_vorticity[index] = old_vorticity[index];
        saved_temperature[index] = old_temperature[index];
        saved_humidity[index] = old_humidity[index];
    }

    if (advection_enabled > 0)
    {
        int32_t mode;
        // 阶段 1：对每个谱模态的垂直列做静力积分。
        for (level = 0; level < level_count; ++level)
        {
            for (mode = 0; mode < spectral_count; ++mode)
            {
                float tendency_sum = 0.0f;
                float temperature_sum = 0.0f;
                int32_t inner;
                for (inner = 0; inner < level_count; ++inner)
                {
                    // G[level, inner]：inner 层温度对目标 level 层位势的贡献。
                    const float coefficient = vertical_g[(size_t)level * (size_t)level_count + (size_t)inner];

                    tendency_sum = tendency_sum + coefficient * explicit_temperature_tendency[field_index(mode, inner, spectral_count)];

                    temperature_sum = temperature_sum + coefficient * saved_temperature[field_index(mode, inner, spectral_count)];
                }

                g_times_tendency[field_index(mode, level, spectral_count)] =
                    tendency_sum;
                g_times_temperature[field_index(mode, level, spectral_count)] =
                    temperature_sum;
            }
        }

        // 阶段 2：逐个水平谱槽求解半隐式 L×L 系统。
        for (level = 0; level < level_count; ++level)
        {
            for (mode = 0; mode < spectral_count; ++mode)
            {
                float result = 0.0f;
                const int32_t degree = wavenumber[mode];
                // n=0 时 λ_n=-n(n+1)=0，不能除以 n(n+1)；常数模态没有重力波，约定 Dbar=0。
                if (degree > 0)
                {
                    const float reciprocal_laplacian =
                        1.0f / (float)(degree * (degree + 1));

                    int32_t inner;
                    for (inner = 0; inner < level_count; ++inner)
                    {
                        const size_t inner_index = field_index(mode, inner, spectral_count);
                        const float reference = reference_temperature[inner];
                        // 构造右端 b_n 的第 inner 个分量。为便于核对公式，把它拆成：
                        //
                        // explicit_tendency = (G R_T'^exp)_inner - T0,inner Abar_exp；
                        // geopotential      = (G T'^{j-1})_inner + T0,inner π^{j-1}；
                        // old_tendency      = R_D^exp/[-λ_n] + Φ_s；
                        // old_divergence    = D^{j-1}/[-λ_n]。
                        //
                        // 随后的 right_hand_side = old_divergence
                        //   + Δτ [ geopotential + old_tendency + Δτ explicit_tendency ]
                        const float explicit_tendency =
                            g_times_tendency[inner_index] -
                            reference * column_pressure_advection[mode];
                        const float geopotential =
                            g_times_temperature[inner_index] +
                            reference * saved_surface_log_pressure[mode];
                        const float old_tendency =
                            explicit_divergence_tendency[inner_index] *
                                reciprocal_laplacian +
                            surface_geopotential[mode];
                        const float old_divergence_value =
                            saved_divergence[inner_index] *
                            reciprocal_laplacian;
                        const float right_hand_side =
                            old_divergence_value + time_step * (geopotential + old_tendency + time_step * explicit_tendency);

                        // inverse_matrix 的布局为 [degree-1][output level][input level]；
                        // 当前取 M_n^{-1}[level, inner]
                        // 使用 degree - 1 是因为 n=0 没有可逆的拉普拉斯项，不构造半隐式矩阵
                        const size_t matrix_index =
                            (size_t)(degree - 1) * (size_t)level_count * (size_t)level_count +
                            (size_t)level * (size_t)level_count +
                            (size_t)inner;

                        result = result + right_hand_side * inverse_matrix[matrix_index];
                    }
                }
                centered_divergence[field_index(mode, level, spectral_count)] = result;
            }
        }
        if (process_id == root_process)
        {
            // 根进程持有 n=0 的实、虚两个槽；它们没有传播重力波，对应 Dbar=0。
            // 即使上面的循环已令 degree==0 时 result 保持零，这里仍显式写零，
            // 既固定该约束，也避免依赖谱槽此前内容或并行分配细节。
            for (level = 0; level < level_count; ++level)
            {
                centered_divergence[field_index(0, level, spectral_count)] =
                    0.0f;
                centered_divergence[field_index(1, level, spectral_count)] =
                    0.0f;
            }
        }

        // 阶段 3：恢复柱质量方程中的完整 Abar。
        // 格点非线性阶段只提供 P 的柱积分 Abar_exp=ΣΔσP；半隐式系统刚求得的
        // Dbar 补上其线性柱积分。因此 Abar 是整柱质量流出率，满足 ∂π/∂τ = -Abar。
        for (mode = 0; mode < spectral_count; ++mode)
        {
            column_mass_divergence[mode] = column_pressure_advection[mode];
        }
        for (level = 0; level < level_count; ++level)
        {
            for (mode = 0; mode < spectral_count; ++mode)
            {
                column_mass_divergence[mode] =
                    column_mass_divergence[mode] +
                    sigma_thickness[level] * centered_divergence[field_index(mode, level, spectral_count)];
            }
        }

        // 阶段 4：恢复完整温度 tendency。
        //
        // 对每个固定的水平谱槽 s=mode，执行同一根 L 层气柱上的矩阵向量积：
        //
        //   R_T',k,s = R_T',k,s^exp - Σ_j Tau_kj · Dbar_j,s。
        //
        // 下标在代码与公式中的对应关系为：
        //   level  <-> k：目标温度层。循环的这一次要写入第 k 层的完整温度 tendency；
        //   mode   <-> s：固定的水平球谐谱槽（某个 (m,n) 系数的实部或虚部）；
        //                    不参与垂直求和，不同 mode 在线性系统中互不混合；
        //   inner  <-> j：垂直求和的来源散度层。第 j 层 Dbar 通过 Tau_kj
        //                    影响目标第 k 层的温度。
        for (level = 0; level < level_count; ++level)
        {
            for (mode = 0; mode < spectral_count; ++mode)
            {
                // 循环结束后：sum = Σ_{j=0}^{L-1} Tau[k,j] · Dbar[j,s]。
                float sum = 0.0f;
                int32_t inner;
                for (inner = 0; inner < level_count; ++inner)
                {
                    sum = sum + vertical_tau[(size_t)level * (size_t)level_count + (size_t)inner] *
                                    centered_divergence[field_index(mode, inner, spectral_count)];
                }
                // complete_temperature_tendency[k,s] = R_T',k,s^exp - sum。
                complete_temperature_tendency[field_index(mode, level, spectral_count)] =
                    explicit_temperature_tendency[field_index(mode, level, spectral_count)] - sum;
            }
        }
    }

    free(g_times_temperature);
    free(g_times_tendency);
}

/*
 * 完成绝热跃蛙更新，并构造 Robert-Asselin 滤波的前两项。
 * short_step_count==0 时先保存
 *   X_filtered = (1-2ν)X^j + νX^{j-1}；
 * 新层 X^{j+1} 在本函数生成，随后 spectrald_finalize 再加 νX^{j+1}，即得到
 *   X~^j = νX^{j-1}+(1-2ν)X^j+νX^{j+1}。
 *
 * column_mass_divergence 保存 Abar=-∂π/∂τ，所以 π 的跃蛙更新带负号。
 * centered_divergence 保存 Dbar，故 D^{j+1}=2Dbar-D^{j-1}；其余变量使用
 * X^{j+1}=X^{j-1}+2Δτ R_X^j。关闭 advection 时则只复制旧层，不推进动力。
 */
void spectrala_finalize(
    int32_t spectral_count, int32_t level_count, int32_t advection_enabled,
    int32_t process_id, int32_t root_process, int32_t short_step_count,
    float double_time_step, float filter_weight, float filter_center_weight,
    const float *column_mass_divergence,
    const float *centered_divergence,
    const float *vorticity_tendency,
    const float *complete_temperature_tendency,
    const float *humidity_tendency,
    const float *saved_surface_log_pressure,
    const float *saved_divergence, const float *saved_vorticity,
    const float *saved_temperature, const float *saved_humidity,
    float *filtered_surface_log_pressure, float *filtered_divergence,
    float *filtered_vorticity, float *filtered_temperature,
    float *filtered_humidity, float *surface_log_pressure,
    float *divergence, float *vorticity, float *temperature, float *humidity)
{
    const int32_t level_values = spectral_count * level_count;
    int32_t index;

    if (short_step_count == 0)
    {
        for (index = 0; index < spectral_count; ++index)
        {
            filtered_surface_log_pressure[index] =
                filter_center_weight * surface_log_pressure[index] +
                filter_weight * saved_surface_log_pressure[index];
        }
        for (index = 0; index < level_values; ++index)
        {
            filtered_divergence[index] =
                filter_center_weight * divergence[index] +
                filter_weight * saved_divergence[index];
            filtered_vorticity[index] =
                filter_center_weight * vorticity[index] +
                filter_weight * saved_vorticity[index];
            filtered_temperature[index] =
                filter_center_weight * temperature[index] +
                filter_weight * saved_temperature[index];
            filtered_humidity[index] =
                filter_center_weight * humidity[index] +
                filter_weight * saved_humidity[index];
        }
    }

    if (advection_enabled > 0)
    {
        for (index = 0; index < spectral_count; ++index)
        {
            surface_log_pressure[index] =
                saved_surface_log_pressure[index] -
                double_time_step * column_mass_divergence[index];
        }
        for (index = 0; index < level_values; ++index)
        {
            divergence[index] =
                2.0f * centered_divergence[index] - saved_divergence[index];

            temperature[index] = double_time_step * complete_temperature_tendency[index] + saved_temperature[index];
            vorticity[index] = double_time_step * vorticity_tendency[index] + saved_vorticity[index];
            humidity[index] = double_time_step * humidity_tendency[index] + saved_humidity[index];
        }
    }
    else
    {
        for (index = 0; index < spectral_count; ++index)
        {
            surface_log_pressure[index] =
                saved_surface_log_pressure[index];
        }
        for (index = 0; index < level_values; ++index)
        {
            divergence[index] = saved_divergence[index];
            temperature[index] = saved_temperature[index];
            vorticity[index] = saved_vorticity[index];
            humidity[index] = saved_humidity[index];
        }
    }

    if (process_id == root_process)
    {
        // 固定全球平均 π 的实、虚槽，避免参考总质量随舍入误差漂移。
        surface_log_pressure[0] = 0.0f;
        surface_log_pressure[1] = 0.0f;
    }
}

/*
 * 把格点物理过程经球谐分析得到的 tendency 显式加到绝热候选新层：
 * X_in=X_*^{j+1}+2Δτ R_phys。函数原地修改状态，输出随后作为谱阻尼的输入。
 */
void spectrald_add_physics(
    int32_t count, float double_time_step,
    const float *vorticity_tendency, const float *temperature_tendency,
    const float *divergence_tendency, const float *humidity_tendency,
    float *vorticity, float *temperature, float *divergence, float *humidity)
{
    int32_t index;
    for (index = 0; index < count; ++index)
    {
        vorticity[index] =
            vorticity[index] + double_time_step * vorticity_tendency[index];
        temperature[index] =
            temperature[index] + double_time_step * temperature_tendency[index];
        divergence[index] =
            divergence[index] + double_time_step * divergence_tendency[index];
        humidity[index] =
            humidity[index] + double_time_step * humidity_tendency[index];
    }
}

/*
 * 根据已经加入物理过程的谱状态计算阻尼 tendency，但暂不更新状态。
 * spectral_diffusion_scale 是 initpm_spectral 构造的高波数选择函数
 * L_{n,k}，不是球谐拉普拉斯特征值 n(n+1)。各阻尼率采用：
 *
 *   R_T = [-(r_T+d_T L)T + r_T T_restore] / [1+2Δτ(r_T+d_T L)]；
 *   R_χ = -d_χ L χ / [1+2Δτ d_χ L]。
 *
 * D 和 η 的瑞利摩擦与谱耗散分别使用各自分母，再把两个 tendency 相加；因此它们
 * 一般不等价于用“摩擦率+耗散率”做一次更新。
 * 谱数组保存绝对涡度 η，根进程对 (m=0,n=1) 模态补回行星涡度，使实际受阻尼的是相对涡度 ζ=η-f。
 */
void spectrald_dissipation(
    int32_t spectral_count, int32_t level_count, int32_t process_id,
    int32_t root_process, int32_t sponge_enabled, float double_time_step,
    float planetary_vorticity, float sponge_rate,
    const float *temperature_drag_rate,
    const float *temperature_diffusion_rate,
    const float *humidity_diffusion_rate, const float *friction_rate,
    const float *divergence_diffusion_rate,
    const float *vorticity_diffusion_rate,
    const float *spectral_diffusion_scale,
    const float *restoration_temperature, const float *divergence,
    const float *vorticity, const float *temperature, const float *humidity,
    float *surface_log_pressure, float *temperature_tendency,
    float *humidity_tendency, float *divergence_tendency,
    float *vorticity_tendency, float *divergence_friction,
    float *divergence_diffusion_tendency, float *vorticity_friction,
    float *vorticity_diffusion_tendency)
{
    int32_t level;
    int32_t mode;

    // 对每层、每个实/虚谱槽计算恢复、扩散和摩擦 tendency。
    for (level = 0; level < level_count; ++level)
    {
        for (mode = 0; mode < spectral_count; ++mode)
        {
            const size_t index = field_index(mode, level, spectral_count);
            const float diffusion_scale = spectral_diffusion_scale[index];

            const float temperature_rate =
                temperature_drag_rate[level] +
                temperature_diffusion_rate[level] * diffusion_scale;

            temperature_tendency[index] =
                ((-temperature_drag_rate[level] - temperature_diffusion_rate[level] * diffusion_scale) * temperature[index] +
                 temperature_drag_rate[level] * restoration_temperature[index]) /
                (1.0f + double_time_step * temperature_rate);

            humidity_tendency[index] =
                -humidity_diffusion_rate[level] * diffusion_scale * humidity[index] /
                (1.0f + double_time_step * humidity_diffusion_rate[level] * diffusion_scale);

            divergence_friction[index] =
                -friction_rate[level] * divergence[index] /
                (1.0f + double_time_step * friction_rate[level]);

            divergence_diffusion_tendency[index] =
                -divergence_diffusion_rate[level] * diffusion_scale * divergence[index] / (1.0f + double_time_step * divergence_diffusion_rate[level] * diffusion_scale);

            divergence_tendency[index] =
                divergence_friction[index] +
                divergence_diffusion_tendency[index];

            vorticity_friction[index] =
                -friction_rate[level] * vorticity[index] /
                (1.0f + double_time_step * friction_rate[level]);

            vorticity_diffusion_tendency[index] =
                -vorticity_diffusion_rate[level] * diffusion_scale * vorticity[index] /
                (1.0f + double_time_step * vorticity_diffusion_rate[level] * diffusion_scale);

            vorticity_tendency[index] =
                vorticity_friction[index] +
                vorticity_diffusion_tendency[index];
        }
    }

    if (process_id == root_process)
    {
        // n=0 的 π 模态固定为零；index=2 是本谱排列中的 (m=0,n=1) 实部。
        surface_log_pressure[0] = 0.0f;
        surface_log_pressure[1] = 0.0f;
        for (level = 0; level < level_count; ++level)
        {
            const size_t index = field_index(2, level, spectral_count);
            const float diffusion_scale = spectral_diffusion_scale[index];
            const float friction_correction =
                planetary_vorticity * friction_rate[level] /
                (1.0f + double_time_step * friction_rate[level]);
            const float diffusion_correction =
                planetary_vorticity * vorticity_diffusion_rate[level] *
                diffusion_scale / (1.0f + double_time_step * vorticity_diffusion_rate[level] * diffusion_scale);

            vorticity_friction[index] =
                vorticity_friction[index] + friction_correction;
            vorticity_diffusion_tendency[index] =
                vorticity_diffusion_tendency[index] + diffusion_correction;
            vorticity_tendency[index] = vorticity_tendency[index] + planetary_vorticity * (friction_rate[level] / (1.0f + double_time_step * friction_rate[level]) + vorticity_diffusion_rate[level] * diffusion_scale / (1.0f + double_time_step * vorticity_diffusion_rate[level] * diffusion_scale));
        }
    }

    if (sponge_enabled > 0)
    {
        // 顶层海绵只作用于温度扰动：显式拉向其下一层，抑制模式顶反射。
        for (mode = 0; mode < spectral_count; ++mode)
        {
            const size_t top = (size_t)mode;
            const size_t below = field_index(mode, 1, spectral_count);
            temperature_tendency[top] =
                sponge_rate * (temperature[below] - temperature[top]) +
                temperature_tendency[top];
        }
    }
}

/*
 * 应用 spectrald_dissipation 产生的阻尼 tendency：X_out=X_in+2Δτ R_damp。
 * 若本次不是启动短步，还把 νX^{j+1} 加入 spectrala_finalize 预存的两项，
 * 至此 Robert-Asselin 滤波的三层组合才完整。π 不受谱阻尼，但仍参与滤波。
 */
void spectrald_finalize(
    int32_t spectral_count, int32_t level_count, int32_t short_step_count,
    float double_time_step, float filter_weight,
    const float *divergence_tendency, const float *vorticity_tendency,
    const float *temperature_tendency, const float *humidity_tendency,
    float *filtered_surface_log_pressure, float *filtered_divergence,
    float *filtered_vorticity, float *filtered_temperature,
    float *filtered_humidity, float *surface_log_pressure,
    float *divergence, float *vorticity, float *temperature, float *humidity)
{
    const int32_t level_values = spectral_count * level_count;
    int32_t index;

    for (index = 0; index < level_values; ++index)
    {
        vorticity[index] =
            vorticity[index] + double_time_step * vorticity_tendency[index];
        temperature[index] =
            temperature[index] + double_time_step * temperature_tendency[index];
        divergence[index] =
            divergence[index] + double_time_step * divergence_tendency[index];
        humidity[index] =
            humidity[index] + double_time_step * humidity_tendency[index];
    }

    if (short_step_count == 0)
    {
        for (index = 0; index < level_values; ++index)
        {
            filtered_vorticity[index] =
                filtered_vorticity[index] + filter_weight * vorticity[index];
            filtered_temperature[index] =
                filtered_temperature[index] + filter_weight * temperature[index];
            filtered_divergence[index] =
                filtered_divergence[index] + filter_weight * divergence[index];
            filtered_humidity[index] =
                filtered_humidity[index] + filter_weight * humidity[index];
        }

        for (index = 0; index < spectral_count; ++index)
        {
            filtered_surface_log_pressure[index] =
                filtered_surface_log_pressure[index] +
                filter_weight * surface_log_pressure[index];
        }
    }
}

/*
 * 构造 σ 交错垂直网格。sigma_half[k] 保存第 k 个整层的下边界
 * σ_{k+1/2}（模式顶 σ=0 隐含，不入数组），最后一个元素应为地表 σ=1。
 *
 * sigma_mode=-1：直接使用配置界面；sigma_mode=1：等厚层；其他值：使用
 * f(s)=3s/4+7s³/4-3s⁴/2 的拉伸网格。随后统一计算
 * Δσ_k、1/(2Δσ_k) 和整层中心 σ_k。调用者负责保证自定义界面严格递增。
 */
void initpm_vertical(
    int32_t level_count, int32_t sigma_mode, const float *configured_sigma,
    float *sigma_half, float *sigma_thickness,
    float *reciprocal_double_sigma, float *sigma_full)
{
    int32_t level;
    if (sigma_mode == -1)
    {
        for (level = 0; level < level_count; ++level)
        {
            sigma_half[level] = configured_sigma[level];
        }
    }
    else if (sigma_mode == 1)
    {
        for (level = 0; level < level_count; ++level)
        {
            sigma_half[level] =
                (float)(level + 1) / (float)level_count;
        }
    }
    else
    {
        for (level = 0; level < level_count; ++level)
        {
            const float sigma =
                (float)(level + 1) / (float)level_count;
            sigma_half[level] = 0.75f * sigma +
                                1.75f * integer_power(sigma, 3) -
                                1.5f * integer_power(sigma, 4);
        }
    }

    // 顶边界为隐含的 0，因此首层厚度就是首个已存下边界。
    sigma_thickness[0] = sigma_half[0];
    for (level = 1; level < level_count; ++level)
    {
        sigma_thickness[level] =
            sigma_half[level] - sigma_half[level - 1];
    }
    for (level = 0; level < level_count; ++level)
    {
        reciprocal_double_sigma[level] =
            0.5f / sigma_thickness[level];
    }
    sigma_full[0] = 0.5f * sigma_half[0];
    for (level = 1; level < level_count; ++level)
    {
        sigma_full[level] =
            0.5f * (sigma_half[level - 1] + sigma_half[level]);
    }
}

static float integer_power(float base, int32_t exponent)
{
    float result = 1.0f;
    int32_t power = exponent;
    float factor = base;
    while (power > 0)
    {
        if ((power & 1) != 0)
        {
            result = result * factor;
        }
        power >>= 1;
        if (power > 0)
        {
            factor = factor * factor;
        }
    }
    return result;
}

/*
 * 初始化逐层阻尼率、谱高波数选择函数和谱索引辅助数组。
 * 所有阻尼时间尺度（restim/tfrc/tdiss*）统一以“天”为单位，代码不做单位
 * 检查，由配置者负责；正的天数 τ 转换为无量纲率 1/(2π τ)，非正值关闭
 * 对应过程。对 degree=n 构造
 *   L_{n,k}=0                                      (n<diffusion_offset)，
 *   L_{n,k}=[(n-diffusion_offset)/(N-diffusion_offset)]^{p_k} (其余)，
 * 因而截断波数 n=N 处 L=1。每个复系数以相邻实、虚两个 float 保存，所以
 * spectral_diffusion_scale、wavenumber 和 spectral_norm 都成对写入。
 */
void initpm_spectral(
    int32_t level_count, int32_t truncation, int32_t diffusion_offset,
    int32_t spectral_count, float two_pi,
    float temperature_scale, const int32_t *diffusion_order,
    const float *restoration_time_days, const float *friction_time_days,
    const float *divergence_diffusion_time_days,
    const float *vorticity_diffusion_time_days,
    const float *temperature_diffusion_time_days,
    const float *humidity_diffusion_time_days,
    const float *reference_temperature_kelvin, float *temperature_drag_rate,
    float *friction_rate, float *divergence_diffusion_rate,
    float *vorticity_diffusion_rate, float *temperature_diffusion_rate,
    float *humidity_diffusion_rate, float *spectral_diffusion_scale,
    int32_t *wavenumber, float *spectral_norm,
    float *nondimensional_reference_temperature)
{
    int32_t level;
    int32_t mode;

    for (level = 0; level < level_count; ++level)
    {
        temperature_drag_rate[level] =
            restoration_time_days[level] > 0.0f
                ? 1.0f / (two_pi * restoration_time_days[level])
                : 0.0f;

        friction_rate[level] =
            friction_time_days[level] > 0.0f
                ? 1.0f / (two_pi * friction_time_days[level])
                : 0.0f;

        divergence_diffusion_rate[level] =
            divergence_diffusion_time_days[level] > 0.0f
                ? 1.0f / (two_pi * divergence_diffusion_time_days[level])
                : 0.0f;

        vorticity_diffusion_rate[level] =
            vorticity_diffusion_time_days[level] > 0.0f
                ? 1.0f / (two_pi * vorticity_diffusion_time_days[level])
                : 0.0f;

        temperature_diffusion_rate[level] =
            temperature_diffusion_time_days[level] > 0.0f
                ? 1.0f / (two_pi * temperature_diffusion_time_days[level])
                : 0.0f;

        humidity_diffusion_rate[level] =
            humidity_diffusion_time_days[level] > 0.0f
                ? 1.0f / (two_pi * humidity_diffusion_time_days[level])
                : 0.0f;

        {
            // normalization 使最高保留波数的选择函数恰为 1。
            const float normalization = 1.0f / integer_power((float)(truncation - diffusion_offset), diffusion_order[level]);

            int32_t zonal;
            mode = 0;
            for (zonal = 0; zonal <= truncation; ++zonal)
            {
                int32_t degree;
                for (degree = zonal; degree <= truncation; ++degree)
                {
                    const float shifted = (float)(degree - diffusion_offset);
                    const float value = degree >= diffusion_offset ? (float)((double)normalization * (double)integer_power(shifted, diffusion_order[level])) : 0.0f;

                    spectral_diffusion_scale[field_index(
                        mode, level, spectral_count)] = value;
                    spectral_diffusion_scale[field_index(
                        mode + 1, level, spectral_count)] = value;
                    mode += 2;
                }
            }
        }
    }

    {
        // 谱排列依次遍历 m=zonal，再遍历 n=degree；每个 (m,n) 占 Re/Im 两槽。
        // norm 随 m 交替符号，保留 PlaSiC 与球谐适配层约定的相位/归一化。
        float norm = 1.0f / sqrtf(2.0f);
        int32_t zonal;
        mode = 0;
        for (zonal = 0; zonal <= truncation; ++zonal)
        {
            int32_t degree;
            for (degree = zonal; degree <= truncation; ++degree)
            {
                wavenumber[mode] = degree;
                wavenumber[mode + 1] = degree;
                spectral_norm[mode] = norm;
                spectral_norm[mode + 1] = norm;
                mode += 2;
            }
            norm = -norm;
        }
    }

    // 有量纲参考温度输入保持为 K；无量纲结果写入单独输出。
    for (level = 0; level < level_count; ++level)
    {
        nondimensional_reference_temperature[level] =
            reference_temperature_kelvin[level] / temperature_scale;
    }
}

/*
 * 为半隐式散度求解预计算每个总波数 n=1..N 的垂直逆矩阵 B_n：
 *
 *   K = G·Tau + T0·Δσ^T,
 *   M_n = Δτ² K + I/[n(n+1)],
 *   B_n = M_n^{-1}。
 *
 * 矩阵只依赖垂直网格、参考态、时间步和 n，
 * 因而初始化时求逆一次，运行时所有同 n 的 m 模态直接复用。
 */
void makebm(
    int32_t level_count, int32_t truncation, float time_step,
    const float *reference_temperature, const float *sigma_thickness,
    const float *vertical_g, const float *vertical_tau,
    float *inverse_matrix)
{
    const float squared_step = time_step * time_step;
    const size_t matrix_values =
        (size_t)level_count * (size_t)level_count;
    float *coupling_matrix =
        (float *)malloc(matrix_values * sizeof(float));
    float *implicit_system_matrix =
        (float *)malloc(matrix_values * sizeof(float));
    int32_t first_level;
    int32_t second_level;
    int32_t degree;
    if (coupling_matrix == NULL || implicit_system_matrix == NULL)
    {
        abort();
    }

    // 先在独立数组中构造与 n 无关的 Δτ²K。
    for (first_level = 0; first_level < level_count; ++first_level)
    {
        for (second_level = 0; second_level < level_count; ++second_level)
        {
            float dot = 0.0f;
            int32_t inner;
            for (inner = 0; inner < level_count; ++inner)
            {
                dot = dot + vertical_g[(size_t)first_level * (size_t)level_count +
                                       (size_t)inner] *
                                vertical_tau[(size_t)inner * (size_t)level_count +
                                             (size_t)second_level];
            }
            coupling_matrix[(size_t)first_level * (size_t)level_count +
                            (size_t)second_level] =
                squared_step *
                (reference_temperature[first_level] *
                     sigma_thickness[second_level] +
                 dot);
        }
    }
    for (degree = 1; degree <= truncation; ++degree)
    {
        float *degree_inverse_matrix =
            inverse_matrix + (size_t)(degree - 1) * matrix_values;
        int32_t level;
        size_t value_index;
        for (value_index = 0; value_index < matrix_values; ++value_index)
        {
            implicit_system_matrix[value_index] = coupling_matrix[value_index];
        }
        for (level = 0; level < level_count; ++level)
        {
            const size_t diagonal = (size_t)level * (size_t)level_count + (size_t)level;
            implicit_system_matrix[diagonal] =
                implicit_system_matrix[diagonal] +
                1.0f / (float)(degree * (degree + 1));
        }

        matrix_inverse(level_count, implicit_system_matrix, degree_inverse_matrix);
    }

    free(implicit_system_matrix);
    free(coupling_matrix);
}

/*
 * 由 σ 网格和参考温度 T0 构造半隐式算法的三个固定垂直算子：
 *
 *   G：静力积分，G·Tv = Φ-Φ_s；
 *   C=W_σ^{-1}G^T W_σ：从模式顶向下的压力速度累积算子；
 *   Tau：散度到线性温度响应的算子，R_T^lin=-Tau·D。
 *
 * G 与 C 的加权转置关系保证离散形式下 压力功和热力能转换使用一致的 Δσ 权重；
 * C 不是 G 的逆。
 * reference_temperature_difference 保存相邻层 T0 的 “下层减上层”，
 * kappa_reference_temperature=κT0。
 */
void initsi(
    int32_t level_count, float kappa,
    const float *reference_temperature, const float *sigma_half,
    const float *sigma_thickness, const float *reciprocal_double_sigma,
    float *kappa_reference_temperature,
    float *reference_temperature_difference, float *vertical_g,
    float *vertical_c, float *vertical_tau)
{
    float *log_ratio = (float *)malloc((size_t)level_count * sizeof(float));
    int32_t level;
    int32_t inner;
    if (log_ratio == NULL)
    {
        abort();
    }

    for (level = 0; level < level_count; ++level)
    {
        kappa_reference_temperature[level] =
            kappa * reference_temperature[level];
    }
    for (level = 0; level < level_count - 1; ++level)
    {
        reference_temperature_difference[level] =
            reference_temperature[level + 1] -
            reference_temperature[level];
    }
    reference_temperature_difference[level_count - 1] = 0.0f;

    // 第 level 层跨越上下界面的对数厚度 ln(σ下/σ上)；顶层用解析极限处理。
    for (level = 1; level < level_count; ++level)
    {
        log_ratio[level] =
            logf(sigma_half[level]) - logf(sigma_half[level - 1]);
    }

    for (level = 0; level < level_count * level_count; ++level)
    {
        vertical_g[level] = 0.0f;
    }

    // 构造上三角静力积分矩阵 G
    vertical_g[0] = 1.0f;
    for (level = 1; level < level_count; ++level)
    {
        vertical_g[(size_t)level * (size_t)level_count + (size_t)level] =
            1.0f - log_ratio[level] * sigma_half[level - 1] /
                       sigma_thickness[level];
        for (inner = 0; inner < level; ++inner)
        {
            vertical_g[(size_t)inner * (size_t)level_count + (size_t)level] =
                log_ratio[level];
        }
    }
    // C_{k j}=(Δσ_j/Δσ_k)G_{j k}，即 W_σ^{-1}G^T W_σ。
    for (level = 0; level < level_count; ++level)
    {
        for (inner = 0; inner < level_count; ++inner)
        {
            vertical_c[(size_t)level * (size_t)level_count + (size_t)inner] =
                vertical_g[(size_t)inner * (size_t)level_count + (size_t)level] *
                (sigma_thickness[inner] / sigma_thickness[level]);
        }
    }

    {
        /*
         * 构造散度—线性温度耦合矩阵 T_D；数组 vertical_tau 就保存该矩阵。
         *
         *     R_T^lin = -T_D D，
         *
         *     (T_D)[k,j] = Γ[k,j] + κ(T0)[k] C[k,j]。
         *
         * 这里 k 是目标温度整层，j 是提供水平散度 D 的整层，Δσ_j 是
         * sigma_thickness[j]。
         *
         * 对内部层 k>=1，Γ 的
         * 离散形式为（最底层省略含 ΔT0[k] 的第二项）：
         *
         *   Γ[k,j] = Δσ_j/(2Δσ_k) * {
         *       ΔT0[k-1] * [σ_(k-1/2) - H(k-1-j)]
         *     + ΔT0[k]   * [σ_(k+1/2) - H(k-j)] }，
         *
         * 其中 ΔT0[k]=T0[k+1]-T0[k]，H(q)=1 (q>=0)，否则为 0。
         *
         * k=0 在下面用其边界极限单独处理。
         */
        float current_delta = reference_temperature_difference[0];
        float current_sigma = sigma_half[0];

        // 顶层（k=0）的 Γ[0,j]，加上仅 j=0 非零的 κT0*C[0,j] 项。
        vertical_tau[0] =
            0.5f * current_delta * (current_sigma - 1.0f) +
            kappa_reference_temperature[0];
        for (inner = 1; inner < level_count; ++inner)
        {
            // 对 j>0，顶层边界极限给出 Γ[0,j]=ΔT0[0]Δσ_j/2。
            vertical_tau[(size_t)inner] =
                0.5f * current_delta * sigma_thickness[inner];
        }
        for (level = 1; level < level_count; ++level)
        {
            // 第 k 层上、下界面附近的参考温度差及对应 σ 界面。
            const float previous_delta = current_delta;
            const float previous_sigma = current_sigma;
            current_delta = reference_temperature_difference[level];
            current_sigma = sigma_half[level];
            for (inner = 0; inner < level_count; ++inner)
            {
                /*
                 * above       = H(k-1-j)：j 层在 k 层上界或其上方；
                 * at_or_above = H(k-j)：  j 层在 k 层内或其上方。
                 *
                 * 这两个量分别用于上、下界面。用 0/1 而非分支写入公式，
                 * 直接实现上述 Γ[k,j] 中的离散阶跃函数。
                 */
                const float at_or_above = inner <= level ? 1.0f : 0.0f;
                const float above = inner < level ? 1.0f : 0.0f;

                // 先累积花括号内的两侧参考温度梯度贡献。
                float value = previous_delta * (previous_sigma - above);
                if (level < level_count - 1)
                {
                    value = value + current_delta * (current_sigma - at_or_above);
                }

                // 乘 Δσ_j/(2Δσ_k)，完成参考温度垂直平流部分 Γ[k,j]。
                value = value * reciprocal_double_sigma[level] * sigma_thickness[inner];

                if (inner <= level)
                {
                    /* C[k,j] 仅在 j<=k 时非零；加上 κT0[k]C[k,j] 的
                     * 绝热压缩/膨胀贡献。这个条件也避免依赖 C 的三角零元。 */
                    value = value + kappa_reference_temperature[level] * vertical_c[(size_t)level * (size_t)level_count + (size_t)inner];
                }
                vertical_tau[(size_t)level * (size_t)level_count + (size_t)inner] = value;
            }
        }
    }

    free(log_ratio);
}

// 由高斯纬度 μ=sinφ 预计算 cos²φ=1-μ² 及 secφ；后者用于 Robert 风还原。
void inilat(
    int32_t latitude_count, const double *sine_latitude,
    float *cosine_squared, float *reciprocal_cosine)
{
    int32_t latitude;
    for (latitude = 0; latitude < latitude_count; ++latitude)
    {
        cosine_squared[latitude] = (float)(1.0 - sine_latitude[latitude] * sine_latitude[latitude]);
        reciprocal_cosine[latitude] =
            1.0f / sqrtf(cosine_squared[latitude]);
    }
}

// 将摩擦回热和谱扩散回热同时累加到现有温度 tendency，不覆盖其他物理源项。
