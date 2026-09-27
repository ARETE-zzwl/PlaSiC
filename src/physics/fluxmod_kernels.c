#include "fluxmod_kernels.h"

#include <math.h>
#include <stddef.h>

/*
 * 将二维场的索引 (horizontal, level) 转换为一维数组索引。
 *
 * 数组内存布局为：
 *
 *     field[level][horizontal]
 *
 * 即同一垂直层的所有水平格点在内存中连续存储。
 *
 * 对应的一维索引为：
 *
 *     index = horizontal + level * nhor
 *
 * 参数：
 *     horizontal ：水平格点编号，范围通常为 [0, nhor-1]
 *     level      ：垂直层编号
 *     nhor       ：水平格点总数
 *
 * 返回：
 *     对应的一维数组下标
 */
static size_t field_index(int32_t horizontal, int32_t level, int32_t nhor)
{
    return (size_t)horizontal + (size_t)level * (size_t)nhor;
}
// int 的位数取决于编译器和平台，虽然现代系统中通常也是 32 位，但标准没有严格保证。
// int32_t 明确保证是 32 位

/*
 * ============================================================================
 * 地表通量交换系数计算
 * flux_surface_exchange_coefficients
 * ============================================================================
 *
 * 由原 Fortran fluxmod.f90 的 surflx 拆分而来；通量驱动顺序现由
 * src/physics/physics_processes.c 的 surface_fluxes() 负责。
 *
 * 该函数计算近地层湍流交换所需的中间变量，包括：
 *
 *     1. 近地面空气温度 dtsa
 *     2. 最低模式层风速平方 zabsu2
 *     3. 最低模式层代表高度 znl
 *     4. 整体 Richardson 数 zri
 *     5. 热量/水汽稳定度修正函数 zrifh
 *     6. 动量稳定度修正函数 zrifm
 *     7. 热量交换系数 dtransh
 *     8. 动量交换系数 dtransm
 */
void flux_surface_exchange_coefficients(
    int32_t nhor, int32_t low_level, int32_t surface_level, int32_t ntsa,
    float vonkarman, float kap, float sigma_low, float rdbrv, float gascon,
    float gravity, float zumin, float vdiff_b, float vdiff_c, float vdiff_d,
    const float *dt, const float *dq, const float *du, const float *dv,
    const float *dz0, const float *dls, float *dtsa, float *dtransh,
    float *dtransm, float *zabsu2, float *znl, float *zri, float *zrifh,
    float *zrifm)
{
    /*
     * 温度从最低 sigma 层外推到地表时使用的指数。
     *
     * kap 通常表示：
     *
     *     κ = R_d / c_p
     *
     * 因此：
     *
     *     zexp = -κ
     */
    const float zexp = -kap;

    /*
     * 最低模式层 sigma 坐标的自然对数：
     *
     *     zlnsig = ln(sigma_low)
     *
     * 由于最低模式层通常满足 0 < sigma_low < 1，
     * 因此 zlnsig 通常小于 0。
     */
    const float zlnsig = logf(sigma_low);

    /*
     * 温度由最低模式层向地表绝热外推的修正因子：
     *
     *     sigma_factor = sigma_low^(-κ)
     */
    const float sigma_factor = powf(sigma_low, zexp);

    /*
     * 虚温中的水汽修正因子。
     *
     *     rdbrv = R_d / R_v
     *
     * 则：
     *
     *     1 / rdbrv - 1 = R_v / R_d - 1
     *
     * 虚温近似可写为：
     *
     *     T_v = T [1 + (R_v/R_d - 1) q]
     */
    const float virtual_factor = 1.0f / rdbrv - 1.0f;

    /*
     * 对每一个水平格点独立计算地表交换系数。
     */
    for (int32_t j = 0; j < nhor; ++j)
    {
        /*
         * low：
         *     当前水平格点 j 在最低模式层上的一维数组下标。
         *
         * surface：
         *     当前水平格点 j 在地表层上的一维数组下标。
         */
        const size_t low = field_index(j, low_level, nhor);
        const size_t surface = field_index(j, surface_level, nhor);

        /*
         * 计算由最低模式层外推得到的近地面空气温度 dtsa。
         *
         * ntsa == 1：
         *     只进行干绝热温度外推：
         *
         *         T_sa = T_low × sigma_low^(-κ)
         *
         * ntsa != 1：
         *     在温度外推的基础上增加虚温水汽修正：
         *
         *         T_v,sa
         *         = T_low × sigma_low^(-κ)
         *         × [1 + (R_v/R_d - 1)q_low]
         */
        if (ntsa == 1)
        {
            dtsa[j] = dt[low] * sigma_factor;
        }
        else
        {
            dtsa[j] = dt[low] * sigma_factor * (1.0f + virtual_factor * dq[low]);
        }

        /*
         * 计算最低模式层水平风速平方：
         *
         *     |V|² = u² + v²
         *
         * 使用 zumin 设置最小值，避免：
         *
         *     1. 后续除以零；
         *     2. 极弱风条件下交换系数完全消失；
         *
         * 注意：
         *     zabsu2 保存的是风速平方，而不是风速本身。
         */
        zabsu2[j] = fmaxf(zumin, du[low] * du[low] + dv[low] * dv[low]);

        /*
         * 估计最低模式层相对于地表的几何高度 znl。
         *
         * 根据静力平衡和测高公式，可以近似写成：
         *
         *     Δz = -(R_d / g) × T_mean × ln(sigma_low)
         *
         * 这里的层平均温度使用：
         *
         *     T_mean = 0.5 × (T_low + T_surface_extrapolated)
         *
         * 因而：
         *
         *     znl = -R_d × 0.5(T_low + T_sa)
         *           × ln(sigma_low) / g
         *
         * 因为 sigma_low < 1，所以 ln(sigma_low) < 0，
         * 最终 znl 通常为正值。
         */
        znl[j] = -gascon * 0.5f * (dt[low] + dtsa[j]) * zlnsig / gravity;

        /*
         * 计算整体 Richardson 数。
         *
         * Richardson 数衡量浮力稳定作用与风切变湍流作用的相对强弱：
         *
         *     Ri ≈ g z ΔT / (T |V|²)
         *
         * Ri > 0：
         *     地表相对较冷，稳定层结，湍流交换受到抑制。
         *
         * Ri < 0：
         *     地表相对较暖，不稳定层结，湍流交换被增强。
         */
        if (ntsa == 1)
        {
            /*
             * 不考虑地表湿度虚温修正：
             *
             *     Ri = g znl (T_sa - T_surface)
             *          / (|V|² T_sa)
             */
            zri[j] = gravity * znl[j] * (dtsa[j] - dt[surface]) / (zabsu2[j] * dtsa[j]);
        }
        else
        {
            /*
             * 考虑地表水汽对虚温的修正：
             *
             *     T_v,surface
             *       = T_surface [1 + (R_v/R_d - 1)q_surface]
             *
             *     Ri = g znl (T_v,sa - T_v,surface)
             *          / (|V|² T_v,sa)
             */
            zri[j] = gravity * znl[j] / (zabsu2[j] * dtsa[j]) * (dtsa[j] - dt[surface] * (1.0f + virtual_factor * dq[surface]));
        }

        /*
         * 计算最低模式层高度与地表粗糙度长度之比：
         *
         *     zbz0 = z / z0
         *
         * dz0[j] 是地表粗糙度长度。
         */
        const float zbz0 = znl[j] / dz0[j];

        /*
         * 中性层结条件下，由对数风廓线得到交换系数的平方根：
         *
         *     neutral = κ_v / ln(1 + z/z0)
         *
         * 其中 κ_v 为 von Kármán 常数。
         *
         * 加 1 的形式可以避免 z/z0 很小时直接计算 ln(0)。
         */
        const float neutral = vonkarman / logf(zbz0 + 1.0f);

        /*
         * 中性条件下的无量纲交换系数：
         *
         *     C_neutral = [κ_v / ln(1 + z/z0)]²
         */
        const float zkblnz2 = neutral * neutral;

        /*
         * 根据 Richardson 数区分不稳定层结和稳定层结。
         */
        if (zri[j] <= 0.0f)
        {
            /*
             * 不稳定或中性层结。
             *
             * 先计算动量与热量稳定度修正公式中的公共分母：
             *
             *     zdenom =
             *       1 + 3 C B sqrt[-Ri(z/z0+1)] C_neutral
             *
             * vdiff_b、vdiff_c 是经验参数。
             */
            const float zdenom = 1.0f + 3.0f * vdiff_c * vdiff_b * sqrtf(-1.0f * zri[j] * (zbz0 + 1.0f)) * zkblnz2;

            /*
             * 不稳定条件下的动量交换稳定度修正：
             *
             *     f_m = 1 - 2 B Ri / zdenom
             *
             * 由于 Ri <= 0，因此第二项通常为正，
             * 即不稳定层结会增强动量交换。
             */
            zrifm[j] = 1.0f - 2.0f * vdiff_b * zri[j] / zdenom;

            /*
             * dls < 1 时采用另一套热量交换参数化。
             *
             * 从使用方式看，该分支用于ocean，
             */
            if (dls[j] < 1.0f)
            {
                /*
                 * 计算地表与近地面空气之间的虚温差：
                 *
                 *     zdth = T_v,surface - T_v,air
                 *
                 * 原表达式写为：
                 *
                 *     -(T_v,air - T_v,surface)
                 */
                const float zdth = -(dtsa[j] - dt[surface] * (1.0f + virtual_factor * dq[surface]));

                /*
                 * 计算浮力对热量交换的增强项：
                 *
                 *     buoyancy =
                 *       0.0016 × zdth^(1/3)
                 *       / |V|
                 *       / C_neutral
                 *
                 * 其中：
                 *
                 *     sqrt(zabsu2) = |V|
                 *
                 * zkblnz2 对应中性交换系数。
                 */
                const float buoyancy = 0.0016f * powf(zdth, 1.0f / 3.0f) / sqrtf(zabsu2[j]) / zkblnz2;

                /*
                 * 特殊地表条件下的热量稳定度修正：
                 *
                 *     f_h = [1 + buoyancy^1.25]^0.8
                 *
                 * 因为：
                 *
                 *     1.25 × 0.8 = 1
                 *
                 * 该形式在浮力项较大时近似呈线性增强。
                 */
                zrifh[j] = powf(1.0f + powf(buoyancy, 1.25f), 0.8f);
            }
            else
            {
                /*
                 * 一般不稳定条件下的热量交换稳定度修正：
                 *
                 *     f_h = 1 - 3 B Ri / zdenom
                 *
                 * 由于 Ri <= 0，因此不稳定层结会增强热量交换。
                 */
                zrifh[j] = 1.0f - 3.0f * vdiff_b * zri[j] / zdenom;
            }
        }
        else
        {
            /*
             * 稳定层结条件 Ri > 0。
             *
             * 稳定层结会抑制垂直湍流交换，因此 zrifh 和 zrifm
             * 通常均小于 1。
             */
            const float zdenom = sqrtf(1.0f + vdiff_d * zri[j]);

            /*
             * 稳定条件下的热量交换修正函数：
             *
             *     f_h =
             *       1 / [1 + 3 B Ri sqrt(1 + D Ri)]
             */
            zrifh[j] = 1.0f / (1.0f + 3.0f * vdiff_b * zri[j] * zdenom);

            /*
             * 稳定条件下的动量交换修正函数：
             *
             *     f_m =
             *       1 / [1 + 2 B Ri / sqrt(1 + D Ri)]
             */
            zrifm[j] = 1.0f / (1.0f + 2.0f * vdiff_b * zri[j] / zdenom);
        }

        /*
         * 计算中性条件下、包含实际风速的交换速度尺度：
         *
         *     ztrans = |V| × C_neutral
         */
        const float ztrans = sqrtf(zabsu2[j]) * zkblnz2;

        /*
         * 加入大气稳定度修正，得到最终的热量和动量交换系数。
         *
         *     dtransh = |V| C_neutral f_h
         *     dtransm = |V| C_neutral f_m
         *
         * 它们通常具有速度量纲，并将传递给后续：
         *
         *     mkshfl  ：感热通量
         *     mkevap  ：蒸发和潜热通量
         *     mkstress：地表动量通量
         */
        dtransh[j] = ztrans * zrifh[j];
        dtransm[j] = ztrans * zrifm[j];
    }
}

/*
 * ============================================================================
 * 地表风应力与摩擦耗散加热
 * flux_mkstress
 * ============================================================================
 *
 * 该函数根据地表动量交换系数 dtransm，计算最低模式层受到的地表摩擦
 * 作用，包括：
 *
 *     1. 隐式更新后的最低模式层纬向风 zun
 *     2. 隐式更新后的最低模式层经向风 zvn
 *     3. 纬向风倾向 dudt
 *     4. 经向风倾向 dvdt
 *     5. 摩擦耗散加热 dtdt
 *     6. 地表纬向风应力 dtaux
 *     7. 地表经向风应力 dtauy
 *     8. 摩擦速度三次方 dust3 = u_*^3
 *
 * 时间离散说明：
 *
 *     deltsec2 表示 2 倍的模式基本积分步长：
 *
 *         deltsec2 = 2 * Δt
 *
 *     这通常对应 leapfrog 时间离散中从旧时间层 t-Δt 更新到
 *     新时间层 t+Δt 的总时间跨度：
 *
 *         (X^{n+1} - X^{n-1}) / (2Δt)
 *
 *     因此，本函数中所有除以 deltsec2 的表达式，表示跨越两个
 *     基本时间步长的平均变化率，而不是除以单个 Δt。
 *
 * 地表摩擦采用半隐式处理。以纬向风为例，其离散形式可理解为：
 *
 *         (u_new - u_old) / (2Δt) = -K u_new
 *
 * 其中：
 *
 *     u_old：
 *         输入数组 du 中当前用于物理过程更新的最低层风速。
 *
 *     u_new：
 *         经过地表摩擦作用后的风速 zun。
 *
 * 整理后得到：
 *
 *         u_new = u_old / (1 + 2Δt K)
 *
 * 在代码中，无量纲系数 zkdiff 已经包含 2Δt，因此写成：
 *
 *         u_new = u_old / (1 + zkdiff)
 *
 * 经向风 v 的计算方式相同。
 *
 * 使用隐式方法的主要原因是：当地表交换较强、近地面风速较大或
 * 模式时间步较长时，显式离散可能导致风速振荡或数值不稳定；
 * 隐式离散能够保证摩擦作用平稳地减小风速。
 */
void flux_mkstress(
    int32_t nhor, int32_t low_level, int32_t surface_level, int32_t ndheat,
    float gravity, float deltsec2, float gascon, float dsigma_low, float cpd,
    float cpv_cpd_minus1, const float *dtransm, const float *dt, const float *dq,
    const float *du, const float *dv, const float *dp, float *dudt,
    float *dvdt, float *dtdt, float *dtaux, float *dtauy, float *dust3,
    float *zun, float *zvn)
{
    /*
     * 计算地表动量交换隐式离散中的公共系数：
     *
     *     zkonst1 = g (2Δt) / (R_d Δσ)
     */
    const float zkonst1 = gravity * deltsec2 / (gascon * dsigma_low);

    /*
     * 将最低模式层风速变化转换为地表动量通量时使用的公共系数：
     *
     *     zkonst2 = Δσ / [(2Δt) g]
     *
     * 乘以地表气压 dp 后：
     *
     *     dp × zkonst2
     *       = p_s Δσ / [(2Δt) g]
     *
     * 它表示最低模式层单位面积空气质量除以跨越两个基本步长的
     * 时间间隔。
     *
     * 注意：
     *
     *     这里的时间间隔是 2Δt，而不是 Δt。
     */
    const float zkonst2 = dsigma_low / deltsec2 / gravity;

    /*
     * 对所有水平格点逐点计算地表风应力。
     *
     * 不同水平格点之间没有数据依赖，因此每个格点可以独立计算。
     */
    for (int32_t j = 0; j < nhor; ++j)
    {
        /*
         * low：
         *     当前水平格点 j 在最低模式层上的一维数组索引。
         *
         * surface：
         *     当前水平格点 j 在地表层上的一维数组索引。
         *
         * 数组采用 level-major 布局：
         *
         *     index = horizontal + level * nhor
         */
        const size_t low = field_index(j, low_level, nhor);
        const size_t surface = field_index(j, surface_level, nhor);

        /*
         * 计算当前格点的无量纲隐式动量交换系数：
         *
         *     zkdiff =
         *       [g (2Δt) / (R_d Δσ)]
         *       × dtransm
         *       / T_surface
         *
         * 即：
         *
         *     zkdiff =
         *       g (2Δt) dtransm
         *       ----------------
         *       R_d Δσ T_surface
         *
         * dtransm 是已经包含风速、中性交换系数和大气稳定度修正的
         * 地表动量交换速度。注意这个方程可以由(3.3)方程推导得到
         *
         * zkdiff 越大，说明：
         *
         *     1. 地表动量交换越强；
         *     2. 地表摩擦对最低层风速的削弱越明显。
         */
        const float zkdiff = zkonst1 * dtransm[j] / dt[surface];

        /*
         * 使用隐式方法更新最低模式层的纬向风：
         *
         *     u_new = u_old / (1 + zkdiff)
         *
         * 其中：
         *
         *     u_old = du[low]
         *     u_new = zun[j]
         *
         * 由于正常情况下 zkdiff >= 0，因此：
         *
         *     |u_new| <= |u_old|
         *
         * 表示地表摩擦削弱纬向风，但不会改变其原始方向。
         */
        zun[j] = du[low] / (1.0f + zkdiff);

        /*
         * 使用相同的隐式方法更新最低模式层的经向风：
         *
         *     v_new = v_old / (1 + zkdiff)
         *
         * 纬向风和经向风使用同一个标量交换系数，因此两个分量按照
         * 相同的比例衰减，水平风向本身不会因这一标量摩擦项而旋转。
         */
        zvn[j] = dv[low] / (1.0f + zkdiff);

        /*
         * 将纬向风变化转换成纬向风倾向：
         *
         *     dudt += (u_new - u_old) / (2Δt)
         *
         * 注意：
         *
         *     deltsec2 = 2Δt
         *
         * 所以这里计算的是从 leapfrog 旧时间层到新时间层之间的
         * 平均纬向风变化率。
         *
         * 使用累加而不是直接赋值，是因为 dudt 中可能已经包含其他
         * 物理过程或动力过程产生的纬向风倾向。
         */
        dudt[low] = dudt[low] + (zun[j] - du[low]) / deltsec2;

        /*
         * 将经向风变化转换成经向风倾向：
         *
         *     dvdt += (v_new - v_old) / (2Δt)
         */
        dvdt[low] = dvdt[low] + (zvn[j] - dv[low]) / deltsec2;

        /*
         * 计算地表摩擦导致的动能耗散，并将耗散的动能转换为热能。
         *
         * 摩擦前的单位质量水平动能为：
         *
         *     KE_old = 0.5 (u_old^2 + v_old^2)
         *
         * 摩擦后的单位质量水平动能为：
         *
         *     KE_new = 0.5 (u_new^2 + v_new^2)
         *
         * 因此，跨越 2Δt 的动能变化率为：
         *
         *     (KE_new - KE_old) / (2Δt)
         *
         * 展开后即为代码中的：
         *
         *     0.5 × [
         *         u_new^2 - u_old^2
         *       + v_new^2 - v_old^2
         *     ] / deltsec2
         *
         * 地表摩擦通常使：
         *
         *     KE_new < KE_old
         *
         * 所以上述动能变化率为负。代码在最前面添加负号，将损失的
         * 动能转换成正的摩擦耗散加热：
         *
         *     heating =
         *       -(KE_new - KE_old)/(2Δt)
         *       --------------------------------
         *       c_pd [1 + cpv_cpd_minus1 q]
         *
         * cpd：
         *     干空气定压比热 c_pd。
         *
         * cpv_cpd_minus1：
         *     水汽对空气热容量的修正系数，通常与
         *
         *         c_pv/c_pd - 1
         *
         *     有关。
         *
         * 1 + cpv_cpd_minus1*dq[low]：
         *     湿空气有效定压比热相对于干空气定压比热的修正因子。
         *
         * 最终 zdtdt 的量纲为温度变化率。
         */
        const float zdtdt = -(zun[j] * zun[j] - du[low] * du[low] + zvn[j] * zvn[j] - dv[low] * dv[low]) / deltsec2 * 0.5f / cpd / (1.0f + cpv_cpd_minus1 * dq[low]);

        /*
         * ndheat 控制是否将地表摩擦耗散的动能反馈为大气加热。
         *
         * ndheat > 0：
         *     将摩擦耗散加热加入最低模式层的温度倾向。
         *
         * ndheat <= 0：
         *     地表摩擦仍然削弱风速，但损失的动能不会反馈到温度场。
         */
        if (ndheat > 0)
        {
            dtdt[low] = dtdt[low] + zdtdt;
        }

        /*
         * 计算地表纬向动量通量，即纬向风应力 dtaux。
         *
         * 代码形式为：
         *
         *     dtaux =
         *       p_s Δσ / [(2Δt)g]
         *       × zkdiff
         *       × u_new
         *
         * 将 zkdiff 的定义代入：
         *
         *     dtaux =
         *       p_s / (R_d T_surface)
         *       × dtransm
         *       × u_new
         *
         * 根据理想气体状态方程：
         *
         *     ρ_surface = p_s / (R_d T_surface)
         *
         * 因此可写成：
         *
         *     dtaux = ρ_surface dtransm u_new
         *
         * 即空气密度、动量交换速度和更新后纬向风速的乘积。
         *
         * dp[j] 从该公式中的使用方式看表示地表气压 p_s。
         */
        dtaux[j] = dp[j] * zkonst2 * zkdiff * zun[j];

        /*
         * 计算地表经向动量通量，即经向风应力 dtauy：
         */
        dtauy[j] = dp[j] * zkonst2 * zkdiff * zvn[j];

        /*
         * 首先计算地表风应力矢量的模：
         *
         *     |τ| = sqrt(dtaux^2 + dtauy^2)
         *
         * 摩擦速度 u_* 定义为：
         *
         *     u_* = sqrt(|τ| / ρ)
         *
         * 利用理想气体状态方程：
         *
         *     1/ρ = R_d T_surface / p_s
         *
         * 得到：
         *
         *     u_*^2 =
         *       |τ| R_d T_surface / p_s
         *
         * 下面这一行首先计算的实际上是 u_*^2，
         * 尽管临时变量名称为 ustar。
         */
        float ustar = sqrtf(dtaux[j] * dtaux[j] + dtauy[j] * dtauy[j]) * gascon * dt[surface] / dp[j];

        /*
         * 对前面得到的 u_*^2 开平方，得到真正的摩擦速度：
         *
         *     u_* = sqrt(u_*^2)
         */
        ustar = sqrtf(ustar);

        /*
         * 保存摩擦速度的三次方：
         *
         *     dust3 = u_*^3
         */
        dust3[j] = ustar * ustar * ustar;
    }
}

/*
 * ============================================================================
 * 地表感热通量与最低模式层温度更新
 * flux_mkshfl
 * ============================================================================
 *
 * 该函数根据地表热量交换系数 dtransh，计算地表与最低模式层大气之间
 * 的感热交换，包括：
 *
 *     1. 隐式更新后的最低模式层温度 ztn
 *     2. 最低模式层温度倾向 dtdt
 *     3. 地表感热通量 dshfl
 *     4. 感热通量对地表温度的导数 dshdt
 *     5. 更新后最低层空气外推到地表的温度 dtsa
 *     6. 最低模式层到地表的温度外推因子 zfac
 *
 * --------------------------------------------------------------------------
 * 时间离散说明
 * --------------------------------------------------------------------------
 *
 * deltsec2 表示两倍的模式基本积分步长：
 *
 *     deltsec2 = 2 * Δt
 *
 * 因此，本函数中除以 deltsec2 的温度变化率应理解为：
 *
 *     温度在 2Δt 时间跨度内的平均变化率
 *
 * --------------------------------------------------------------------------
 * 最低模式层温度向地表外推
 * --------------------------------------------------------------------------
 *
 * 模式中的 dt[low] 是最低 sigma 层上的温度，而地表感热交换需要使用
 * 与地表处于同一参考高度的空气温度。
 *
 * 代码使用干绝热关系，将最低模式层温度外推到地表：
 *
 *     T_air,surface = T_low * sigma_low^(-κ)
 *
 * 其中：
 *
 *     κ = R_d / c_p
 *
 * --------------------------------------------------------------------------
 * 感热交换的隐式离散
 * --------------------------------------------------------------------------
 *
 * 地表感热通量通常可以写为：
 *
 *     H = ρ c_p C_h |V| (T_air,surface - T_surface)
 *
 * 本函数中的 dtransh 已经综合包含：
 *
 *     1. 近地面风速；
 *     2. 中性条件下的交换系数；
 *     3. 大气稳定度修正。
 *
 * 因此可将感热交换简写为：
 *
 *     H ∝ ρ c_p dtransh
 *         (zfac * T_low - T_surface)
 *
 * 为提高数值稳定性，通量中的最低层温度采用新时间层温度 ztn，
 * 而不是旧时间层温度 dt[low]。其隐式离散形式可理解为：
 *
 *     T_new - T_old
 *       = zkdiff * (T_surface - zfac * T_new)
 *
 * 整理后得到：
 *
 *                     T_old + zkdiff * T_surface
 *     T_new = ------------------------------------------------
 *                  1 + zkdiff * zfac
 *
 * 这正是代码中 ztn 的计算公式。
 */
void flux_mkshfl(
    int32_t nhor, int32_t low_level, int32_t surface_level, float kap,
    float sigma_low, float gravity, float deltsec2, float gascon,
    float dsigma_low, float cpd, float cpv_cpd_minus1, const float *dtransh,
    const float *dt, const float *dq, const float *dp, float *dtdt,
    float *dshfl, float *dshdt, float *dtsa, float *ztn, float *zfac)
{
    /*
     * 温度从最低 sigma 层向地表进行干绝热外推时使用的指数：
     *
     *     zexp = -κ
     *
     * 其中：
     *
     *     κ = R_d / c_p
     *
     * kap 即通常所说的 Poisson 常数。
     */
    const float zexp = -kap;

    /*
     * 计算感热交换隐式离散中的公共系数：
     *
     *     zkonst1 = g (2Δt) / (R_d Δσ)
     */
    const float zkonst1 = gravity * deltsec2 / (gascon * dsigma_low);

    /*
     * 将最低模式层温度变化转换为地表能量通量时使用的公共系数：
     *
     *     zkonst2 = Δσ / [(2Δt) g]
     */
    const float zkonst2 = dsigma_low / deltsec2 / gravity;

    /*
     * 计算最低模式层温度向地表进行干绝热外推的修正因子：
     *
     *     extrapolation = sigma_low^(-κ)
     *
     * 由于最低模式层通常满足：
     *
     *     0 < sigma_low < 1
     *
     * 且 κ > 0，因此该修正因子通常略大于 1。
     *
     * 这表示在干绝热条件下，从较高的最低模式层向较低的地表外推时，
     * 空气温度通常会升高。
     */
    const float extrapolation = powf(sigma_low, zexp);

    /*
     * 对每一个水平格点独立计算地表感热交换。
     */
    for (int32_t j = 0; j < nhor; ++j)
    {
        /*
         * low：
         *     当前水平格点 j 在最低模式层上的一维数组索引。
         *
         * surface：
         *     当前水平格点 j 在地表层上的一维数组索引。
         *
         * 数组采用如下线性存储方式：
         *
         *     index = horizontal + level * nhor
         */
        const size_t low = field_index(j, low_level, nhor);
        const size_t surface = field_index(j, surface_level, nhor);

        /*
         * 计算当前格点的无量纲隐式热交换系数：
         *
         *     zkdiff =
         *       [g (2Δt) / (R_d Δσ)]
         *       * dtransh
         *       / T_surface
         *
         * 即：
         *
         *                 g (2Δt) dtransh
         *     zkdiff = -------------------------
         *                R_d Δσ T_surface
         *
         */
        const float zkdiff = zkonst1 * dtransh[j] / dt[surface];

        /*
         * 计算水汽对空气定压比热的修正因子：
         *
         *     moisture_heat_capacity = max(1, 1 + cpv_cpd_minus1 * q)
         *
         * 湿空气的定压比热可近似表示为：
         *
         *     c_p,moist = c_pd * (1 + cpv_cpd_minus1 * q)
         *
         * 其中：
         *
         *     cpd：
         *         干空气定压比热 c_pd。
         *
         *     cpv_cpd_minus1：
         *         水汽对定压比热的修正系数，通常与
         *
         *             c_pv / c_pd - 1
         *
         *         有关。
         *
         *     dq[low]：
         *         最低模式层比湿。
         *
         * 使用 fmaxf 将修正因子限制为不小于 1，避免由于异常比湿值
         * 或数值误差造成湿空气有效热容量小于干空气热容量。
         */
        const float moisture_heat_capacity = fmaxf(1.0f, 1.0f + cpv_cpd_minus1 * dq[low]);

        /*
         * 保存最低模式层温度向地表进行干绝热外推的因子：
         *
         *     zfac = sigma_low^(-κ)
         *
         * 该因子本身与水平格点无关，但代码仍逐格点写入 zfac，
         * 便于后续诊断或与原 Fortran 数组接口保持一致。
         */
        zfac[j] = extrapolation;

        /*
         * 隐式求解经过地表感热交换后的最低模式层温度。
         *
         * 记：
         *
         *     T_old     = dt[low]
         *     T_new     = ztn[j]
         *     T_surface = dt[surface]
         *     F         = zfac[j]
         *
         * 隐式离散关系可写为：
         *
         *     T_new - T_old
         *       = zkdiff * (T_surface - F * T_new)
         *
         * 右侧使用的是新时间层最低层温度外推到地表后的值：
         *
         *     T_air,surface,new = F * T_new
         *
         * 展开：
         *
         *     T_new - T_old
         *       = zkdiff * T_surface
         *         - zkdiff * F * T_new
         *
         * 将包含 T_new 的项移到左侧：
         *
         *     T_new * (1 + zkdiff * F)
         *       = T_old + zkdiff * T_surface
         *
         * 最终得到：
         *
         *                 T_old + zkdiff * T_surface
         *     T_new = ----------------------------------
         *                    1 + zkdiff * F
         *
         * 对应代码：
         *
         *     ztn =
         *       (dt[low] + zkdiff * dt[surface])
         *       / (1 + zkdiff * zfac)
         *
         * 这种写法保证强感热交换条件下，最低层温度能够稳定地趋近于
         * 与地表温度相一致的状态。
         */
        ztn[j] = (dt[low] + zkdiff * dt[surface]) / (1.0f + zkdiff * zfac[j]);

        /*
         * 将感热交换造成的最低层温度变化转换为温度倾向：
         *
         *     dtdt += (T_new - T_old) / (2Δt)
         *
         * 因为：
         *
         *     deltsec2 = 2Δt
         *
         * 所以这里不是除以单个模式基本步长，而是除以 leapfrog
         * 更新跨越的完整时间跨度。
         */
        dtdt[low] = dtdt[low] + (ztn[j] - dt[low]) / deltsec2;

        /*
         * 计算地表感热通量 dshfl。
         *
         * 代码中的表达式为：
         *
         *     dshfl =
         *       [Δσ / ((2Δt)g)]
         *       * zkdiff
         *       * p_s
         *       * c_pd
         *       * moisture_heat_capacity
         *       * (T_air,surface,new - T_surface)
         *
         * 其中：
         *
         *     T_air,surface,new = ztn * zfac
         *
         * 将 zkonst2 和 zkdiff 展开：
         *
         *     zkonst2 * zkdiff * p_s
         *
         *       = [Δσ / ((2Δt)g)]
         *         * [g(2Δt)dtransh / (R_d Δσ T_surface)]
         *         * p_s
         *
         *       = p_s dtransh / (R_d T_surface)
         *
         * 根据理想气体状态方程：
         *
         *     ρ_surface = p_s / (R_d T_surface)
         *
         * 因此感热通量可改写为：
         *
         *     dshfl =
         *       ρ_surface
         *       * dtransh
         *       * c_p,moist
         *       * (T_air,surface,new - T_surface)
         *
         * 其中：
         *
         *     c_p,moist =
         *       c_pd * moisture_heat_capacity
         *
         * 这就是常见的整体空气动力学感热通量公式。
         *
         * 虽然原始代码表达式中显式出现 deltsec2 = 2Δt，
         * 但将 zkdiff 展开后，2Δt 会完全约去。因此最终计算出的
         * 物理感热通量本身不直接依赖模式时间步长。
         *
         * 通量正负号取决于：
         *
         *     ztn * zfac - dt[surface]
         */
        dshfl[j] = zkonst2 * zkdiff * dp[j] * cpd * moisture_heat_capacity * (ztn[j] * zfac[j] - dt[surface]);

        /*
         * 计算感热通量对地表温度的线性化导数 dshdt。
         *
         * 代码采用：
         *
         *     dshdt =
         *       -zkonst2
         *       * zkdiff
         *       * p_s
         *       * c_pd
         *       * moisture_heat_capacity
         *
         * 将公共系数展开后可写为：
         *
         *     dshdt =
         *       -ρ_surface
         *       * dtransh
         *       * c_p,moist
         *
         * 负号来自感热通量中的温差项：
         *
         *     T_air,surface - T_surface
         *
         * 在保持其他量不变时：
         *
         *     ∂(T_air,surface - T_surface)
         *     -------------------------------- = -1
         *             ∂T_surface
         *
         * dshdt 通常用于地表能量平衡方程的隐式求解或线性化迭代。
         */
        dshdt[j] = -1.0f * zkonst2 * zkdiff * dp[j] * cpd * moisture_heat_capacity;

        /*
         * 使用隐式更新后的最低模式层温度，重新计算外推到地表高度的
         * 近地面空气温度：
         *
         *     dtsa = T_new * sigma_low^(-κ)
         *
         * 即：
         *
         *     dtsa = ztn * zfac
         *
         * 后续地表能量平衡或其他地表通量计算可以使用该更新后的
         * 近地面空气温度。
         */
        dtsa[j] = ztn[j] * zfac[j];
    }
}

/*
 * ============================================================================
 * 地表蒸发、最低模式层比湿更新与潜热通量
 * flux_mkevap
 * ============================================================================
 *
 * 该函数根据地表热量/水汽交换系数 dtransh，计算地表与最低模式层
 * 大气之间的水汽交换，包括：
 *
 *     1. 隐式更新后的最低模式层比湿 zqn
 *     2. 最低模式层比湿倾向 dqdt
 *     3. 地表蒸发水质量通量 devap
 *     4. 地表潜热通量 dlhfl
 *     5. 潜热通量对地表温度的线性化导数 dlhdt
 *     6. 根据地表可用水量限制实际蒸发量
 *
 * 最低模式层比湿倾向采用：
 *
 *     (q_new - q_old) / (2Δt)
 *
 * 因此代码中的 dqdt 使用 deltsec2。
 *
 * 地表可用水量 dwatc 是在一个基本时间步长 Δt 内被消耗或释放的
 * 水储量，因此代码中水量限制使用 deltsec，而不是 deltsec2。
 *
 * 这一点非常重要：
 *
 *     deltsec2 用于大气状态变量的 leapfrog 时间倾向；
 *
 *     deltsec  用于地表水储量在单个模式时间步内的消耗限制。
 *
 * --------------------------------------------------------------------------
 * 水汽交换的整体空气动力学形式
 * --------------------------------------------------------------------------
 *
 * 地表水汽通量通常可以写为：
 *
 *     E = ρ C_h |V| wetfac (q_air - q_surface)
 *
 * 其中：
 *
 *     ρ：
 *         近地面空气密度。
 *
 *     C_h |V|：
 *         地表热量和水汽交换速度。
 *
 *     wetfac：
 *         土壤湿润度影响的蒸发效率因子。
 *
 *     q_surface：
 *         地表饱和比湿或地表边界比湿。
 *
 *     q_air：
 *         最低模式层空气比湿。
 *
 * 本函数中的 dtransh 已经综合包含：
 *
 *     1. 近地面风速；
 *     2. 地表粗糙度；
 *     3. 大气稳定度修正。
 *
 * dwetfac 则进一步限制地表实际能够提供水汽的能力。
 *
 * --------------------------------------------------------------------------
 * 最低模式层比湿的隐式更新
 * --------------------------------------------------------------------------
 *
 * 隐式离散关系可写为：
 *
 *     q_new - q_old
 *       = zkdiff (q_surface - q_new)
 *
 * 整理后：
 *
 *     q_new (1 + zkdiff)
 *       = q_old + zkdiff q_surface
 *
 * 因此：
 *
 *                 q_old + zkdiff q_surface
 *     q_new = --------------------------------
 *                       1 + zkdiff
 *
 * 这正是代码中 implicit_q 的计算公式。
 *
 * 随后代码使用：
 *
 *     q_new = max(q_old, implicit_q)
 *
 * 强制该物理过程只能增加最低层比湿，而不能降低最低层比湿。
 * 因此该函数处理的是地表蒸发，不处理露水沉降或向下水汽通量。
 */
void flux_mkevap(
    int32_t nhor, int32_t low_level, int32_t surface_level, float gravity,
    float deltsec, float deltsec2, float gascon, float dsigma_low, float tmelt,
    float lv, float ls, float ra2, float ra4, const float *dtransh,
    const float *dt, const float *dq, const float *dp, const float *dwetfac,
    const float *dls, const float *dwatc, float *dqdt, float *devap,
    float *dlhfl, float *dlhdt, float *zqn)
{
    /*
     * 计算水汽交换隐式离散中的公共系数：
     *
     *     zkonst1 = g (2Δt) / (R_d Δσ)
     */
    const float zkonst1 = gravity * deltsec2 / (gascon * dsigma_low);

    /*
     *     zkonst2 = Δσ / [(2Δt) g]
     */
    const float zkonst2 = dsigma_low / deltsec2 / gravity;

    /*
     * 对每一个水平格点独立计算地表蒸发和潜热交换。
     */
    for (int32_t j = 0; j < nhor; ++j)
    {
        /*
         * low：
         *     当前水平格点 j 在最低模式层上的一维数组索引。
         *
         * surface：
         *     当前水平格点 j 在地表层上的一维数组索引。
         *
         * 数组线性索引形式为：
         *
         *     index = horizontal + level * nhor
         */
        const size_t low = field_index(j, low_level, nhor);
        const size_t surface = field_index(j, surface_level, nhor);

        /*
         * 计算当前格点的无量纲隐式水汽交换系数：
         *
         *     zkdiff =
         *       dwetfac
         *       * [g (2Δt) / (R_d Δσ)]
         *       * dtransh
         *       / T_surface
         *
         * 即：
         *
         *                 dwetfac g (2Δt) dtransh
         *     zkdiff = --------------------------------
         *                  R_d Δσ T_surface
         *
         * 其中：
         *
         *     dwetfac[j]：
         *         dwetfac 越小，说明地表水分供应越受限，蒸发越弱。
         *
         *     dtransh[j]：
         *         热量和水汽交换速度，已经包含风速、粗糙度和
         *         大气稳定度修正。
         *
         *     dt[surface]：
         *         地表温度。
         */
        const float zkdiff = dwetfac[j] * zkonst1 * dtransh[j] / dt[surface];

        /*
         * 隐式求解经过地表水汽交换后的最低模式层比湿。
         *
         * 记：
         *
         *     q_old     = dq[low]
         *     q_new     = implicit_q
         *     q_surface = dq[surface]
         *
         * 隐式关系为：
         *
         *     q_new - q_old
         *       = zkdiff (q_surface - q_new)
         *
         * 展开：
         *
         *                 q_old + zkdiff q_surface
         *     q_new = --------------------------------
         *                       1 + zkdiff
         *
         * 当 zkdiff 较大时，q_new 会更接近地表边界比湿 q_surface；
         * 当 zkdiff 较小时，q_new 更接近原始空气比湿 q_old。
         */
        const float implicit_q = (dq[low] + zkdiff * dq[surface]) / (1.0f + zkdiff);

        /*
         * 强制更新后的最低模式层比湿不小于原始比湿：
         *
         *     zqn = max(q_old, implicit_q)
         */
        zqn[j] = fmaxf(dq[low], implicit_q);

        /*
         * 当 dls[j] > 0 时，根据地表实际可用水量 dwatc 对蒸发量
         * 施加上限。
         *
         * 从代码逻辑看，dls 用于区分不同地表类型。
         *
         * 对于需要考虑有限地表水储量的格点，蒸发不能超过当前时间步
         * 内地表能够提供的水量。
         */
        if (dls[j] > 0.0f)
        {
            /*
             * 根据地表可用水量计算最低模式层比湿允许达到的最大值。
             *
             * 首先：
             *
             *     dwatc[j] / deltsec
             *
             * 表示在一个模式基本时间步 Δt 内，地表最多能够提供的
             * 单位时间水量。
             *
             * 注意这里使用 deltsec = Δt，而不是 deltsec2 = 2Δt。
             *
             * 这是因为 dwatc 是地表水储量约束，其消耗发生在一个
             * 实际模式时间步内；而 dqdt 的大气状态更新时间离散采用
             * leapfrog 的 2Δt 跨度。
             *
             * 后面的单位换算：
             *
             *     * 1000 / dp[j] / zkonst2
             *
             * 将地表可用水质量通量转换成最低模式层允许增加的比湿。
             
             * dwatc / deltsec：
             *     一个基本模式步长内允许的最大蒸发水深通量，单位 m s^-1。
             * 乘以 1000：
             *     利用水密度约 1000 kg m^-3，将水深通量转换为
             *     水质量通量 kg m^-2 s^-1。
             * dp * zkonst2：
             *     最低模式层单位面积空气质量除以 2Δt，
             *     单位同样为 kg m^-2 s^-1。
             * 两者相除得到无量纲的最大允许比湿增量 Δq_max。
             
             * 最大新比湿 available_q。
             * 因此：
             *     available_q
             *       = q_old + 地表水储量允许的最大比湿增量
             */
            const float available_q = dwatc[j] / deltsec * 1000.0f / dp[j] / zkonst2 + dq[low];

            /*
             * 实际更新后的最低模式层比湿不能超过地表可用水量所允许的
             * 最大值：
             *
             *     q_new = min(q_new, available_q)
             */
            zqn[j] = fminf(zqn[j], available_q);
        }

        /*
         * 将最低模式层比湿变化转换成比湿倾向：
         *
         *     dqdt += (q_new - q_old) / (2Δt)
         */
        dqdt[low] = dqdt[low] + (zqn[j] - dq[low]) / deltsec2;

        /*
         * 根据最低模式层比湿变化，反推出地表蒸发水质量通量。
         *
         * 代码形式为：
         *
         *     devap =
         *       -p_s Δσ / [(2Δt)g]
         *       * (q_new - q_old)
         *       / 1000
         *
         * 即：
         *
         *     devap =
         *       -dp * zkonst2
         *       * Δq
         *       / 1000
         *
         * 其中：
         *
         *     p_s Δσ/g：
         *         最低模式层单位面积空气质量。
         *
         *     q_new - q_old：
         *         由地表蒸发引起的空气比湿增加。
         *
         * 因此：
         *
         *     p_s Δσ/g * (q_new-q_old)
         *
         * 表示最低模式层单位面积水汽质量的增加量。
         *
         * 再除以 2Δt，得到相应的水质量通量。
         *
         * 前面的负号来自 PlaSiC 对蒸发通量的符号约定：
         *
         *     当 q_new > q_old 时，devap < 0。
         *
         * 因此在本代码约定中，地表向大气蒸发对应负的 devap。
         *
         */
        devap[j] = -dp[j] * zkonst2 / 1000.0f * (zqn[j] - dq[low]);

        /*
         * 根据地表温度和地表类型，选择水汽相变潜热。
         *
         * 判断条件为：
         *
         *     T_surface > T_melt
         *
         * 或：
         *
         *     dls < 0.5
         *
         * 满足任一条件时：
         *
         *     latent_heat = lv
         *
         * 即使用液态水的汽化潜热。
         *
         * 否则：
         *
         *     latent_heat = ls
         *
         * 即使用冰雪的升华潜热。
         *
         * 参数含义：
         *
         *     tmelt：
         *         水的融点温度。
         *
         *     lv：
         *         latent heat of vaporization，汽化潜热。
         *
         *     ls：
         *         latent heat of sublimation，升华潜热。
         *
         *     dls：
         *  陆海掩膜, 1 = 陆地, 0 = 海洋
         */
        const float latent_heat = (dt[surface] > tmelt || dls[j] < 0.5f)
                                      ? lv
                                      : ls;

        /*
         * 将蒸发水质量通量转换成地表潜热通量：
         *
         *     dlhfl = devap * latent_heat * 1000
         *
         * 即：
         *
         *     LH = E L
         *
         * 其中：
         *
         *     E：
         *         蒸发水质量通量。
         *
         *     L：
         *         汽化潜热或升华潜热。
         *
         * 乘以 1000 与 devap 计算中除以 1000 的单位换算相对应。
         *
         * 由于蒸发时 devap 通常为负，因此 dlhfl 的符号也通常为负。
         * 这表示地表通过蒸发损失能量
         */
        dlhfl[j] = devap[j] * latent_heat * 1000.0f;

        /*
         * 计算饱和比湿经验公式中使用的温度差：
         *
         *     temperature_delta = T_surface - ra4
         *
         * ra4 是饱和水汽压或饱和比湿经验公式中的常数。
         *
         * 后续潜热通量对地表温度的导数中包含：
         *
         *     1 / (T_surface - ra4)^2
         */
        const float temperature_delta = dt[surface] - ra4;

        /*
         * 计算潜热通量对地表温度的线性化导数 dlhdt。
         *
         * 地表饱和比湿通常是地表温度的强非线性函数。
         * 这里使用的经验关系，其温度导数具有如下结构：
         *
         *     dq_surface/dT_surface
         *       =
         *       ra2 (tmelt - ra4) q_surface
         *       --------------------------------
         *          (T_surface - ra4)^2
         *
         * 因此潜热通量对地表温度的近似导数为：
         *
         *     dlhdt =
         *       -L
         *       * zkdiff
         *       * p_s Δσ / [(2Δt)g]
         *       * dq_surface/dT_surface
         *
         * 对应代码：
         *
         *     dlhdt =
         *       -latent_heat
         *       * zkdiff
         *       * zkonst2
         *       * dp
         *       * ra2
         *       * (tmelt-ra4)
         *       * q_surface
         *       / (T_surface-ra4)^2
         *
         * 将 zkonst2 和 zkdiff 展开后：
         *
         *     zkonst2 * zkdiff * dp
         *
         * 中的 2Δt 会相消。因此 dlhdt 作为通量对温度的导数，
         * 最终并不直接依赖模式时间步长。
         */
        dlhdt[j] = -1.0f * latent_heat * zkdiff * zkonst2 * dp[j] * ra2 * (tmelt - ra4) * dq[surface] / (temperature_delta * temperature_delta);

        /*
         * 如果实际潜热通量 dlhfl 恰好为零，则将潜热通量对地表温度的
         * 导数也强制设为零：
         *
         *     if (dlhfl == 0)
         *         dlhdt = 0
         *
         * 这样可以保证：
         *
         *     1. 没有实际蒸发时，不向地表能量方程加入潜热反馈；
         *
         *     2. 水量限制或湿度限制使蒸发停止时，线性化导数与实际
         *        通量状态保持一致；
         *
         *     3. 避免出现“潜热通量为零，但其温度导数非零”的数值
         *        不一致。
         *
         * 注意这里使用浮点数精确相等判断：
         *
         *     dlhfl[j] == 0.0f
         *
         * 这是原始实现的一部分，此处保持不变。
         */
        if (dlhfl[j] == 0.0f)
        {
            dlhdt[j] = 0.0f;
        }
    }
}
