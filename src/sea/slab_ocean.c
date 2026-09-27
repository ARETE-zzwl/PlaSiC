#include "slab_ocean.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

/*
 * 初始化单层混合层海洋的参数和时间步长。
 *
 * 参数说明：
 *   nhor    : 水平网格点数。
 *   ntspd   : 每个太阳日内的混合时间步数。
 *   solar_day: 一个太阳日对应的模型时间长度。
 *   taunc   : 输入/输出；输入通常是以“天”为尺度的时间常数，输出转换为模型时间单位。Newtonian cooling timescale 
 *   layer_depth : 混合层层厚，对所有水平网格点相同。
 *   ymld    : 输出；将 layer_depth 扩展到所有水平网格点。
 *   dtmix   : 输出；一次混合时间步对应的模型时间长度。
 */
void ocean_initialize_parameters(
    int32_t nhor, int32_t ntspd, float solar_day,
    float *taunc, float layer_depth, float *ymld, float *dtmix)
{
    /*
     * 单层混合层：层厚 layer_depth 对所有水平点相同，
     * 因此直接复制成混合层厚度场 ymld(horizontal)。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        ymld[horizontal] = layer_depth;
    }

    /* 一个太阳日均分为 ntspd 个混合时间步。 */
    *dtmix = solar_day / (float)ntspd;

    /* 将以“太阳日”为单位给出的时间常数转换为模型时间单位。 */
    const float cooling_timescale_in_days = *taunc;
    const float cooling_timescale_in_model_time =
        solar_day * cooling_timescale_in_days;
    *taunc = cooling_timescale_in_model_time;
}

/*
 * 在两个相邻月的水平场之间做线性插值。
 *
 * 参数说明：
 *   month1/month2 : 两个时间切片的零基索引。
 *   weight        : 第二个时间切片的权重；第一个切片权重为 1 - weight。
 *   monthly       : 按 [month][horizontal] 连续存储的输入场。
 *   field         : 输出的当前时刻水平场。
 *
 * 当 weight = 0 时完全取 month1；weight = 1 时完全取 month2。
 * 函数本身不限制 weight 范围，因此调用端也可以用它进行线性外推。
 */
void ocean_interpolate_cycle(
    int32_t nhor, int32_t month1, int32_t month2, float weight,
    const float *monthly, float *field)
{
    /* 两个月份在连续数组中的起始偏移。 */
    const size_t first = (size_t)month1 * (size_t)nhor;
    const size_t second = (size_t)month2 * (size_t)nhor;

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* 标准线性插值：(1-w) * x1 + w * x2。 */
        field[horizontal] =
            (1.0f - weight) * monthly[first + (size_t)horizontal] + weight * monthly[second + (size_t)horizontal];
    }
}

/*
 * 将海洋内部状态复制到对外输出数组。
 *
 * 参数说明：
 *   ysst/ymld: 内部海表温度和混合层厚度场，均为单层水平场。
 *   yiflux     : 内部冰相关热通量。
 *   psst/pmld/piflux: 对外输出。
 *
 * 注意：piflux 对 yiflux 取负，说明内部与外部接口采用相反的通量正方向约定。
 */
void ocean_copy_outputs(
    int32_t nhor, const float *ysst, const float *ymld,
    const float *yiflux, float *psst, float *pmld, float *piflux)
{
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        psst[horizontal] = ysst[horizontal];
        pmld[horizontal] = ymld[horizontal];

        /* 转换通量符号约定。 */
        piflux[horizontal] = -yiflux[horizontal];
    }
}

/*
 * 将耦合接口输入复制到海洋模块内部工作数组，并构造冰上总雪/冰水当量。
 *
 * yicesnow 的表达式把冰厚乘以冰密度后除以 1000，再与雪量相加；
 * 若 crhoi 以 kg/m^3、厚度以 m 表示，则该项对应以水密度 1000 kg/m^3 归一化的水当量。
 */
void ocean_copy_inputs(
    int32_t nhor, float crhoi, const float *picec, const float *piced,
    const float *pheat, const float *ppme, const float *ptaux,
    const float *ptauy, const float *pust3,
    const float *psnow, const float *pcliced, float *yicec, float *yiced,
    float *yheat, float *ypme, float *ytaux, float *ytauy, float *yust3,
    float *yicesnow, float *ycliced)
{
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* 以下为逐水平点的一一复制，不做缩放或筛选。 */
        yicec[horizontal] = picec[horizontal];
        yiced[horizontal] = piced[horizontal];
        yheat[horizontal] = pheat[horizontal];
        ypme[horizontal] = ppme[horizontal];
        ytaux[horizontal] = ptaux[horizontal];
        ytauy[horizontal] = ptauy[horizontal];
        yust3[horizontal] = pust3[horizontal];

        /* 冰的水当量与积雪量相加，得到冰雪总量相关内部场。 */
        yicesnow[horizontal] = piced[horizontal] * crhoi / 1000.0f + psnow[horizontal];

        ycliced[horizontal] = pcliced[horizontal];
    }
}

/*
 * 参数说明：
 *   dtmix   : 当前混合时间步长度。
 *   crhos   : 海水密度。
 *   cps     : 海水比热容。
 *   tfreeze : 冻结温度。
 *   yls     : 海陆掩膜，yls < 1 表示参与海洋计算。
 *   yiced   : 海冰厚度或冰量；<= 0 视为无冰，> 0 视为有冰。
 *   yheat   : 表面热通量项。
 *   yfldo   : 另一个进入海洋表层热收支的通量项。 fl：flux，通量 do：deep ocean，深海
 *   ymld    : 混合层厚度场；ymld[horizontal] 为该点的混合层厚度。
 *   ysst    : 单精度输入温度场。
 *   zsst    : 双精度工作温度场。
 *   yiflux  : 冰相关热通量工作数组。
 */
void ocean_mksst_begin(
    int32_t nhor, float dtmix, float crhos, float cps,
    float tfreeze, const float *yls, const float *yiced,
    const float *yheat, const float *yfldo, const float *ymld,
    const float *ysst, double *zsst, float *yiflux)
{
    /*
     * 热通量转换为温度变化时的公共系数：dt / (rho * cp)。
     * 后续还会除以混合层厚度 ymld，从而得到温度增量。
     */
    const float temperature_response_factor = dtmix / (crhos * cps);

    /* 处理每个水平点的表层热收支和初始冰通量。 */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* 混合层温度从 float 复制到 double，后续计算采用双精度工作变量。 */
        zsst[horizontal] = (double)ysst[horizontal];

        /* 每个时间步开始时先清零冰通量工作量。 */
        yiflux[horizontal] = 0.0f;

        /*
         * 海点且无冰：净热通量 yheat + yfldo 直接改变表层海温。
         * 温度增量 = 通量 * dt / (rho * cp * 层厚)。
         */
        if (yls[horizontal] < 1.0f && yiced[horizontal] <= 0.0f)
        {
            const float increment =
                (yheat[horizontal] + yfldo[horizontal]) *
                temperature_response_factor / ymld[horizontal];
            zsst[horizontal] += (double)increment;
        }

        /*
         * 海点且有冰：不在此处直接改变海温，而把两项热通量的相反数记录到 yiflux。
         */
        if (yls[horizontal] < 1.0f && yiced[horizontal] > 0.0f)
        {
            yiflux[horizontal] = -yfldo[horizontal] - yheat[horizontal];
        }
    }

    /*
     * 保留原始 MKIFLUX 的条件语义：条件检查的是输入数组 ysst，而不是已更新的 zsst。
     * 这一区别可能影响刚刚因热通量跨越冻结点的网格，review 时不可擅自替换。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        if (ysst[horizontal] < tfreeze && yls[horizontal] < 1.0f)
        {
            /*
             * 将 zsst 调整到 tfreeze 所需的能量换算为通量并累加到 yiflux：
             *   ΔF = (tfreeze - zsst) * ymld / temperature_response_factor。
             * 随后把表层工作温度强制设置为冻结点。
             */
            yiflux[horizontal] =
                (float)((double)yiflux[horizontal] +
                        ((double)tfreeze - zsst[horizontal]) *
                            (double)ymld[horizontal] /
                            (double)temperature_response_factor);
            zsst[horizontal] = (double)tfreeze;
        }
    }
}

/*
 * 完成海表温度（SST）时间步的收尾处理。
 *
 * 参数说明：
 *   nocean  - 海洋模式开关；大于 0 表示使用预报海洋。
 *   nfluko  - 海洋热通量修正方案编号；0 表示关闭，1 表示在线
 *              计算 Newtonian restoring 修正通量。
 *   dtmix   - 本次混合/海洋时间步长度。
 *   crhos   - 海水密度。
 *   cps     - 海水定压比热容。
 *   tfreeze - 海水冻结温度。
 *   yls     - 海陆掩膜；代码以小于 1.0 表示海洋格点。
 *   yiced   - 海冰厚度或等价冰量；大于 0 表示存在海冰。
 *   ymld    - 各格点的混合层厚度。
 *   zsst    - 双精度 SST 工作数组。
 *   yiflux  - 海冰相关热通量，函数内进行累加。
 *   ysst    - 单精度 SST 状态数组。
 */
void ocean_mksst_finish(
    int32_t nhor, int32_t nocean, int32_t nfluko,
    float dtmix, float crhos, float cps, float tfreeze,
    const float *yls, const float *yiced, const float *ymld,
    double *zsst, float *yiflux, float *ysst)
{
    /*
     * 温度响应系数：时间步 /（密度 × 比热）。
     */
    const float temperature_response_factor = dtmix / (crhos * cps);
    /*
     * 只有在预报海洋开启且不使用海洋热通量修正时，才在此处显式处理
     * 有海冰格点的冻结点约束。
     */
    if (nocean > 0 && nfluko == 0)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            /* 仅处理“存在海冰”的海洋格点，陆地格点保持不变。 */
            if (yiced[horizontal] > 0.0f && yls[horizontal] < 1.0f)
            {
                /*
                 * 反推出把当前表层温度 zsst 调整到 tfreeze 所需的热通量。
                 * ymld/temperature_response_factor 将温差换算回单位面积通量。
                 */
                const float flux_to_freezing_temperature =
                    (float)(((double)tfreeze - zsst[horizontal]) *
                            (double)ymld[horizontal] /
                            (double)temperature_response_factor);
                /*
                 * 再把该通量换算成温度增量。该写法与上面的反推公式配对，
                 * 理想算术下会使表层温度到达冻结点。
                 */
                const float increment =
                    flux_to_freezing_temperature *
                    temperature_response_factor / ymld[horizontal];
                zsst[horizontal] += (double)increment;
                /* 将本次冻结点调整对应的通量并入海冰热通量诊断。 */
                yiflux[horizontal] += flux_to_freezing_temperature;
            }
        }
    }
    /*
     * zsst 用双精度完成中间计算；时间步结束时写回单精度 ysst。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        ysst[horizontal] = (float)zsst[horizontal];
    }
}

/*
 * 应用 Newtonian SST restoring 通量，并在海洋显热与海冰相变热之间
 * 重新分配能量。
 *
 * 对每个水平格点，本函数使用同一个混合层热收支关系：
 *
 *   C_A = rho_w * c_p * h_mld
 *   dT  = F * dtmix / C_A
 *
 * 其中 C_A 是单位面积混合层热容量，F 是进入海洋混合层的热通量。
 * 因而，通量与温度变化之间可以互相换算：
 *
 *   temperature increment = F * dtmix / (crhos * cps * ymld)
 *   flux needed for dT    = dT * crhos * cps * ymld / dtmix
 *
 * `yclsst2` 的名字中虽然带有数字 2，但它不是第二层或数组下标 2；
 * 它表示“当前模型时刻插值后的气候态 SST”。原始的逐月气候态表是
 * `yclsst`，`yclsst2` 是从该表得到的当前时刻二维场。
 *
 * 根据海冰状态，本函数有三条路径：
 *
 *   1. 无海冰：
 *      restoring 通量全部改变 SST。
 *
 *   2. 有气候态海冰（ycliced > 0）：
 *      将 yclsst2 作为局部目标温度，只施加不会越过该温度的部分；
 *      未用于改变 SST 的通量通过 yifluxr 转交给海冰模块。
 *
 *   3. 有实际海冰但没有气候态海冰：
 *      先将 restoring 通量全部施加到 SST；若 SST 高于海水冻结点，
 *      再计算融化现有海冰所需的潜热，并从海洋中扣除。
 *
 * 第二遍循环还会消耗已有的正 yiflux，将高于约束温度的 SST 冷却下来；
 * 最后强制海洋 SST 不低于 tfreeze，并把对应的冻结/结冰能量记入 yifluxr （r 表示 residual，即“残余量”。）。
 * yiflux 使用海洋侧符号，之后传给海冰模块时会在 ocean_copy_outputs() 中反号。
 */
static void add_flux_correction(
    int32_t nhor, float dtmix, float taunc, float crhos, float crhoi, float cps,
    float clfi, float tfreeze, const float *yls, const float *yclsst2,
    const float *ymld, const float *yiced, const float *ycliced,
    const float *yheat, float *ysst, float *yiflux, float *yifluxr)
{
    /*
     * 热通量到温度变化的公共系数：
     *   dtmix / (crhos * cps)
     * 对单位面积通量 F 还需再除以混合层厚度 ymld，才能得到 dT。
     */
    const float temperature_response_factor = dtmix / (crhos * cps);
    /* 体积热容量 rho_w * c_p，用于从温差反推单位面积热通量。 */
    const float volumetric_heat_capacity = crhos * cps;
    /*
     * zsst：表层 SST 的双精度工作副本。一个时间步内可能连续执行
     *       restoring、融冰冷却和冻结点截断，因此使用双精度减少累计舍入误差。
     *
     * constraint_temperature：第二遍循环使用的局部温度上限/目标：
     *       - 默认是 tfreeze；
     *       - 气候态海冰格点改为当前气候态 SST yclsst2。
     */
    double *zsst = malloc((size_t)nhor * sizeof(*zsst));
    float *constraint_temperature =
        malloc((size_t)nhor * sizeof(*constraint_temperature));

    /* 内存分配失败时清理并终止，避免继续使用空指针。 */
    if (zsst == NULL || constraint_temperature == NULL)
    {
        free(zsst);
        free(constraint_temperature);
        abort();
    }

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* yifluxr 只记录本次修正产生的残余通量，不继承上一个时间步。 */
        yifluxr[horizontal] = 0.0f;
        constraint_temperature[horizontal] = tfreeze;
        zsst[horizontal] = (double)ysst[horizontal];

        /* 陆地格点不参与海洋热通量修正。 */
        if (yls[horizontal] >= 1.0f)
        {
            continue;
        }

        /*
         * 根据当前 SST 偏差即时计算 restoring 通量。
         *
         * taunc > 0 时：
         *
         *   F_R = rho_w*c_p*h*(T_clim - T_sst) / tau
         *
         * 这是方程 dT/dt = (T_clim - T_sst)/tau 乘以混合层单位面积
         * 热容量后的结果。正通量使海洋升温，负通量使海洋降温。
         */
        const float flux_correction =
            taunc > 0.0f ?
                volumetric_heat_capacity *
                    (yclsst2[horizontal] - ysst[horizontal]) *
                    ymld[horizontal] / taunc :
                volumetric_heat_capacity / dtmix *
                    (yclsst2[horizontal] - ysst[horizontal]) *
                    ymld[horizontal] -
                    yheat[horizontal];

        /*
         * 无海冰：restoring 通量全部进入海洋显热。
         *
         *   T <- T + F_R * dtmix / (rho_w*c_p*h)
         */
        if (yiced[horizontal] <= 0.0f)
        {
            const float increment =
                flux_correction * temperature_response_factor /
                ymld[horizontal];
            zsst[horizontal] += (double)increment;
        }
        /*
         * 气候态海冰：yclsst2 是当前时刻的气候态 SST，作为局部目标温度。
         * 只有“朝向目标温度”的 restoring 通量才允许改变 SST，并且不能越过
         * 目标温度；剩余通量转为海冰侧的残余热通量。
         */
        else if (ycliced[horizontal] > 0.0f)
        {
            constraint_temperature[horizontal] = yclsst2[horizontal];
            /*
             * 当前 SST 高于气候态温度，且 restoring 为冷却方向。
             * 到达目标温度所需的通量为：
             *
             *   F_target = rho_w*c_p*h*(yclsst2 - T) / dtmix < 0
             *
             * fmax(F_target, F_R) 使实际冷却量不超过 restoring 提供的量，
             * 从而不会把 SST 冷却到 yclsst2 以下。
             */
            if (zsst[horizontal] > (double)yclsst2[horizontal] && flux_correction < 0.0f)
            {
                const float flux_to_climatology =
                    (float)(((double)yclsst2[horizontal] -
                             zsst[horizontal]) *
                            (double)ymld[horizontal] /
                            (double)temperature_response_factor);
                /* 实际用于改变 SST 的通量：取两者中较大的负值。 */
                const float applied_flux =
                    fmaxf(flux_to_climatology, flux_correction);
                const float increment =
                    applied_flux * temperature_response_factor /
                    ymld[horizontal];
                zsst[horizontal] += (double)increment;
                /*
                 * 未用于 SST 的 restoring 部分转交海冰：
                 *
                 *   yifluxr = -(F_R - F_applied)
                 *
                 * 这里的负号来自海洋侧 yiflux 与海冰侧通量的符号约定。
                 */
                yifluxr[horizontal] =
                    yifluxr[horizontal] -
                    (flux_correction - applied_flux);
            }
            /*
             * 当前 SST 低于气候态温度，且 restoring 为增暖方向。
             * 此时 F_target > 0，使用 fmin(F_target, F_R)，避免升温越过目标。
             */
            else if (zsst[horizontal] < (double)yclsst2[horizontal] && flux_correction > 0.0f)
            {
                const float flux_to_climatology =
                    (float)(((double)yclsst2[horizontal] -
                             zsst[horizontal]) *
                            (double)ymld[horizontal] /
                            (double)temperature_response_factor);
                /* 实际只使用到达目标温度所需的那一部分正通量。 */
                const float applied_flux =
                    fminf(flux_to_climatology, flux_correction);
                const float increment =
                    applied_flux * temperature_response_factor /
                    ymld[horizontal];
                zsst[horizontal] += (double)increment;
                yifluxr[horizontal] =
                    yifluxr[horizontal] -
                    (flux_correction - applied_flux);
            }
            /*
             * 其他方向组合会把 SST 推离目标温度：
             * 不让 restoring 改变 SST，而把全部通量交给海冰侧处理。
             */
            else
            {
                yifluxr[horizontal] -= flux_correction;
            }
        }
        /*
         * 实际海冰但无气候态海冰：先按普通混合层热收支施加 restoring。
         * 与气候态海冰分支不同，这里不把 yclsst2 当作 SST 的硬目标。
         */
        else
        {
            const float increment =
                flux_correction * temperature_response_factor /
                ymld[horizontal];
            zsst[horizontal] += (double)increment;
            /*
             * 若施加后 SST 高于冻结点，则将当前冰量对应的融化潜热
             * 从海洋显热中扣除：
             *
             *   F_melt = rho_i*h_i*L_i / dtmix
             *
             * 这是单位面积、单位时间的等效融冰通量。它转入 yifluxr，
             * 后续传给海冰模块后会改变海冰潜热储量。
             */
            if (zsst[horizontal] > (double)tfreeze)
            {
                /*
                 * 单位面积融冰能量约为冰量 × 冰密度 × 融化潜热；
                 * 除以 dtmix 后得到等效热通量。
                 */
                const float melt_flux =
                    yiced[horizontal] * crhoi * clfi / dtmix;
                const float cooling =
                    melt_flux * temperature_response_factor /
                    ymld[horizontal];
                zsst[horizontal] -= (double)cooling;
                yifluxr[horizontal] -= melt_flux;
            }
        }
    }

    /*
     * 第二遍：处理之前已经存在的 yiflux，并执行最终冻结点约束。
     * 第一遍主要决定 restoring 通量如何分配；第二遍则把已有海冰通量
     * 与 SST 状态重新协调，避免同一份能量同时被计入海洋和海冰。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /*
         * 若 SST 高于局部约束温度，且 yiflux > 0，尝试用已有海冰通量冷却海洋。
         * 到达约束温度所需的通量为：
         *
         *   F_target = rho_w*c_p*h*(T_constraint - T) / dtmix < 0.
         *
         * 由于 yiflux 是可用的正通量，实际抵消量不能超过 -yiflux，
         * 所以使用 max(F_target, -yiflux)。
         */
        if (zsst[horizontal] > (double)constraint_temperature[horizontal] &&
            yiflux[horizontal] > 0.0f && yls[horizontal] < 1.0f)
        {
            const float flux_to_constraint_temperature =
                (float)(((double)constraint_temperature[horizontal] -
                         zsst[horizontal]) *
                        (double)ymld[horizontal] /
                        (double)temperature_response_factor);
            /* 限制实际冷却量，不能消耗超过已有的 yiflux。 */
            const float applied_flux =
                fmaxf(flux_to_constraint_temperature,
                      -yiflux[horizontal]);
            const float increment =
                applied_flux * temperature_response_factor /
                ymld[horizontal];
            zsst[horizontal] += (double)increment;
            yiflux[horizontal] += applied_flux;
        }
        /*
         * 最终硬性冻结点约束：SST 不能低于 tfreeze。
         * 若 T < tfreeze，需要补充：
         *
         *   F_freeze = rho_w*c_p*h*(tfreeze - T) / dtmix > 0.
         *
         * 这部分能量使海水回到冻结点，并作为结冰/冻结残余通量转交海冰。
         */
        if (zsst[horizontal] < (double)tfreeze && yls[horizontal] < 1.0f)
        {
            yifluxr[horizontal] =
                (float)((double)yifluxr[horizontal] +
                        ((double)tfreeze - zsst[horizontal]) *
                            (double)ymld[horizontal] /
                            (double)temperature_response_factor);
            zsst[horizontal] = (double)tfreeze;
        }
        /* 将双精度工作值写回单精度表层 SST。 */
        ysst[horizontal] = (float)zsst[horizontal];
    }

    /* 释放本函数申请的临时数组。 */
    free(zsst);
    free(constraint_temperature);
}

/*
 * 完成一个海洋时间步：生成/应用热通量修正、累积诊断量，并按输出周期
 * 计算时间平均。
 *
 * 返回值：
 *   0 - 当前步不是输出步，诊断量继续累积；
 *   1 - 当前步到达输出周期，诊断累积量已除以样本数得到平均值。
 *
 * 关键控制参数：
 *   nocean = 0：使用规定的气候态 SST，不预报海洋温度；
 *   nocean > 0：使用预报海洋。
 *
 *   nfluko = 0：不使用海洋热通量修正；
 *   nfluko = 1：根据气候态 SST 与当前 SST 的偏差在线计算
 *               Newtonian restoring 修正通量。
 *
 * 诊断累积数组：
 *   yheata  - 表面热通量累计；
 *   yifluxa - 海冰热通量累计；
 *   yfldoa  - 海洋附加/动力热通量累计。
 *
 * naccuout 记录自上次重置以来累计的样本数。
 */
int32_t ocean_finish_step(
    int32_t nhor, int32_t nocean, int32_t nfluko, int32_t nstep,
    int32_t nout, float dtmix, float taunc, float crhos, float crhoi,
    float cps, float clfi, float tfreeze, const float *yls,
    const float *yclsst2, const float *ymld, const float *yiced,
    const float *ycliced, const float *yheat, const float *yfldo,
    float *ysst, float *yiflux, float *yifluxr,
    float *yheata, float *yifluxa, float *yfldoa, int32_t *naccuout)
{
    /*
     * 无预报海洋：直接把 SST 设为气候态值。
     */
    if (nocean == 0)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            if (yls[horizontal] < 1.0f)
            {
                ysst[horizontal] = yclsst2[horizontal];
            }
        }
    }
    /* 预报海洋且启用方案 1：在线计算并应用恢复通量。 */
    else if (nfluko == 1)
    {
        add_flux_correction(
            nhor, dtmix, taunc, crhos, crhoi, cps, clfi, tfreeze, yls,
            yclsst2, ymld, yiced, ycliced, yheat, ysst, yiflux, yifluxr);
    }

    /*
     * 启用在线修正方案时，将 add_flux_correction 产生的残余通量
     * 合并到主海冰热通量 yiflux。
     */
    if (nocean > 0 && nfluko == 1)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            if (yls[horizontal] < 1.0f)
            {
                yiflux[horizontal] += yifluxr[horizontal];
            }
        }
    }

    /*
     * 累积本时间步的各类诊断量。这里不区分海陆，保持调用方提供数组的
     * 原始约定；陆地值应由上游逻辑保证为合适的值。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        yheata[horizontal] += yheat[horizontal];
        yfldoa[horizontal] += yfldo[horizontal];
        yifluxa[horizontal] += yiflux[horizontal];
    }
    /* 本次时间步已贡献一个诊断样本。 */
    *naccuout += 1;

    /*
     * 只有当 nstep 能被 nout 整除时才生成输出平均值。
     * 调用方应保证 nout 非零。
     */
    if (nstep % nout != 0)
    {
        return 0;
    }
    /* 用累计样本数把各诊断累计和转换为时间平均值。 */
    const float count = (float)*naccuout;
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        yheata[horizontal] /= count;
        yfldoa[horizontal] /= count;
        yifluxa[horizontal] /= count;
    }
    /* 返回 1 通知调用方：本步已形成可输出的平均诊断量。 */
    return 1;
}

/*
 * 重置海洋诊断累计量，为下一个输出累计窗口做准备。
 *
 * 所有逐格点累计数组清零，同时把累计样本计数 naccuout 置为 0。
 * 通常应在 finish_step 返回 1、诊断量完成输出之后调用。
 */
void ocean_reset_diagnostics(
    int32_t nhor, float *yheata, float *yifluxa,
    float *yfldoa, int32_t *naccuout)
{
    /* 逐水平格点清空全部诊断累计数组。 */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        yheata[horizontal] = 0.0f;
        yifluxa[horizontal] = 0.0f;
        yfldoa[horizontal] = 0.0f;
    }
    /* 清空累计样本数，使下一输出窗口从零开始计数。 */
    *naccuout = 0;
}
