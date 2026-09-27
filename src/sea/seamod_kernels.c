#include "seamod_kernels.h"

#include <math.h>
#include <stddef.h>

static size_t field_index(int32_t horizontal, int32_t level, int32_t nhor)
{
    return (size_t)horizontal + (size_t)level * (size_t)nhor;
}

/*
 * 初始化海洋和海冰相关状态。
 *
 * 该函数遍历所有水平格点，仅处理满足 dls[j] < 0.5 的海洋格点。
 * 对海洋格点执行以下操作：
 *
 *   1. 将海洋模块提供的状态 c* 复制到大气侧状态 d*。
 *   2. 将海表温度写入大气温度场的指定表面层。
 *   3. 保存表面层当前比湿。
 *   4. 对于非重启动初始化，重新计算：
 *        - 表面饱和比湿
 *        - 相对湿度控制参数
 *        - 表面反照率
 *        - 表面粗糙度
 *
 * 参数命名约定：
 *   c* - 通常表示来自耦合器、海洋或海冰模块的输入状态
 *   d* - 通常表示大气模式侧使用的诊断量或状态量
 *
 * 注意：
 *   dls[j] < 0.5 被用作海洋格点判据。
 *
 * 主要参数说明：
 *   nhor         - 水平格点总数；
 *   surface_level- 大气三维场中的表面层编号；
 *   nrestart     - 重启动标志，等于 1 时保留重启动场中的湿度、反照率等量；
 *   rdbrv        - 干空气与水汽气体常数之比 R_d / R_v；
 *   ra1、ra2、ra4、tmelt、psurf
 *                - 用于计算表面饱和湿度的热力学参数；
 *   albsea、albice
 *                - 开阔海面和海冰的反照率参数；
 *   dz0sea、dz0ice
 *                - 开阔海面和海冰的表面粗糙度；
 *   dwetfacsea、dwetfacice
 *                - 开阔海面和海冰的相对湿度控制参数；
 *   c*           - 海洋/海冰模块传入的状态量；
 *   d*           - 大气模式侧需要初始化或更新的状态量。
 */
void sea_initialize(
    int32_t nhor, int32_t surface_level, int32_t nrestart,
    float rdbrv, float ra1, float ra2, float ra4, float tmelt, float psurf,
    float albsea, float albice, float dz0sea, float dz0ice,
    float dwetfacsea, float dwetfacice,
    const float *dls, const float *cts, const float *csst,
    const float *cmld, const float *cicec, const float *ciced,
    const float *csnow, float *dts, float *dicec, float *diced,
    float *dsnow, float *dsst, float *dmld, float *dt, float *dq,
    float *dqs, float *dwetfac, float *dalb, float *dz0)
{
    /*
     * 湿度转换修正系数。
     *
     * 饱和湿度公式先得到混合比 r，利用
     *
     *     q = r / [1 - (R_v / R_d - 1) r]
     *
     * 将其转换为比湿 q。由于 rdbrv = R_d / R_v，
     * 这里的 R_v / R_d - 1 写作 1 / rdbrv - 1。
     */
    const float humidity_correction = 1.0f / rdbrv - 1.0f;

    /* 遍历所有水平格点。 */
    for (int32_t j = 0; j < nhor; ++j)
    {
        /*
         * 仅处理海洋格点。
         */
        if (dls[j] < 0.5f)
        {
            /*
             * 获取当前水平格点在大气三维场表面层中的一维索引。
             */
            const size_t surface = field_index(j, surface_level, nhor);

            /*
             * 从耦合输入状态初始化大气侧海洋/海冰状态。
             */
            dts[j] = cts[j];     /* 表面温度 */
            dicec[j] = cicec[j]; /* 海冰覆盖率 */
            diced[j] = ciced[j]; /* 海冰厚度 */
            dsnow[j] = csnow[j]; /* 海冰上的积雪状态 */
            dsst[j] = csst[j];   /* 海表温度 */
            dmld[j] = cmld[j];   /* 海洋混合层深度 */

            /*
             * 用海洋表面温度更新大气温度场的表面层。
             */
            dt[surface] = dts[j];

            /*
             * 暂时读取并保存大气表面层原有比湿。
             *
             * 在重启动情况下，后续不会重新计算比湿，因此保留重启动场中的值。
             */
            dqs[j] = dq[surface];

            /*
             * 非重启动初始化时，根据当前表面温度重新计算表面湿度、
             * 相对湿度参数、反照率和粗糙度。
             * nrestart == 1 时保留重启动文件中已有的相关状态。
             */
            if (nrestart != 1)
            {
                /*
                 * 根据表面温度计算饱和湿度。
                 */
                const float saturation_mixing_ratio =
                    rdbrv * ra1 *
                    expf(ra2 * (dt[surface] - tmelt) /
                         (dt[surface] - ra4)) /
                    psurf;
                const float saturation_specific_humidity =
                    saturation_mixing_ratio /
                    (1.0f - humidity_correction *
                                saturation_mixing_ratio);

                /*
                 * 同时更新二维表面饱和比湿和三维大气表面层比湿。
                 */
                dqs[j] = saturation_specific_humidity;
                dq[surface] = saturation_specific_humidity;

                /*
                 * 根据海冰覆盖率对开阔海面和海冰表面的湿度参数做线性混合：
                 *
                 *   dicec = 0：完全采用海面值 dwetfacsea
                 *   dicec = 1：完全采用海冰值 dwetfacice
                 */
                dwetfac[j] = dwetfacsea * (1.0f - dicec[j]) + dwetfacice * dicec[j];

                /*
                 * 根据海冰覆盖率计算格点平均反照率。
                 *
                 * 开阔海面部分使用固定海洋反照率 albsea。
                 *
                 * 海冰部分使用：
                 *
                 *     min(albice, 0.5 + 0.025 * (273 - Tsurface))
                 *
                 * 即海冰反照率随表面温度降低而增大，但不超过 albice。
                 */
                dalb[j] = albsea * (1.0f - dicec[j]) + dicec[j] * fminf(albice, 0.5f + 0.025f * (273.0f - dts[j]));

                /*
                 * 根据海冰覆盖率，对海面粗糙度和海冰粗糙度做线性混合。
                 */
                dz0[j] = dz0sea * (1.0f - dicec[j]) + dz0ice * dicec[j];
            }
        }
    }
}

/*
 * 累积大气向海洋/海冰耦合模块提供的通量，并在耦合时刻计算平均值。
 *
 * 每次调用时，对海洋格点累积以下量：
 *
 *   cheata - 总热通量
 *   cpmea  - 降水、蒸发相关的淡水通量 （pme：precipitation minus evaporation，即“降水减蒸发”）
 *   cprsa  - 固态降水
 *   ctauxa - 纬向或第一方向表面应力
 *   ctauya - 经向或第二方向表面应力
 *   cust3a - 风速或相关三次方量
 *   cshfla - 感热通量
 *   cshdta - 感热通量对温度的导数
 *   clhfla - 潜热通量
 *   clhdta - 潜热通量对温度的导数
 *   cswfla - 短波辐射通量
 *   clwfla - 长波辐射通量
 *
 * naccua 记录当前累计了多少个大气时间步。
 *
 * 当同时满足以下条件时，进入耦合时刻：
 *
 *   1. nstep 能被 ncpl_atmos_ice 整除；
 *   2. nkits == 0。
 *
 * 到达耦合时刻后，所有累计量除以 naccua，转换成时间平均值。
 *
 * 参数说明：
 *   ncpl_atmos_ice- 大气与海洋/海冰模块之间的耦合间隔（时间步数）；
 *   nkits         - 初始显式 Euler 启动步数；为 0 时表示进入正常积分阶段
 *
 * 累计数组在本函数中不会被清零；它们应在上一个耦合周期结束时由
 * sea_finish_coupling 清零。
 *
 * 返回值：
 *   0 - 当前不是实际耦合时刻，输出数组中仍保存累计和
 *   1 - 当前是耦合时刻，输出数组已经转换为时间平均值
 */
int32_t sea_prepare_step(
    int32_t nhor, int32_t surface_level, int32_t nstep,
    int32_t ncpl_atmos_ice, int32_t nkits, int32_t *naccua,
    const float *dls, const float *dshfl, const float *dswfl,
    const float *dlwfl, const float *dlhfl, const float *dprl,
    const float *dprc, const float *devap, const float *dprs,
    const float *dtaux, const float *dtauy, const float *dust3,
    const float *dshdt, const float *dlhdt, float *cheata, float *cpmea,
    float *cprsa, float *ctauxa, float *ctauya, float *cust3a, float *cshfla,
    float *cshdta, float *clhfla, float *clhdta, float *cswfla,
    float *clwfla)
{
    /*
     * 将当前大气时间步产生的通量加入累计数组。
     */
    for (int32_t j = 0; j < nhor; ++j)
    {
        /* 仅对海洋格点累积海气耦合量。 */
        if (dls[j] < 0.5f)
        {
            /*
             * 短波和长波辐射通量存储在带垂直维度的场中，
             * 因此需要取指定表面层。
             */
            const size_t surface = field_index(j, surface_level, nhor);

            /*
             * 累积总表面热通量。
             * 当前时间步的热通量由感热、短波、长波和潜热四部分组成。
             * dswfl 和 dlwfl 是三维场，因此取当前水平格点的表面层值；
             * dshfl 和 dlhfl 已经是二维表面量，可以直接使用。
             */
            cheata[j] = dshfl[j] + dswfl[surface] + dlwfl[surface] + dlhfl[j] + cheata[j];

            /*
             * 累积淡水通量：液态降水、对流降水和蒸发共同构成水量收支项。
             * 各分量的正负号约定沿用大气通量计算模块，此处只负责求和。
             */
            cpmea[j] = dprl[j] + dprc[j] + devap[j] + cpmea[j];

            /* 累积固态降水。 */
            cprsa[j] = dprs[j] + cprsa[j];

            /* 累积两个水平方向的表面应力。 */
            ctauxa[j] = dtaux[j] + ctauxa[j];
            ctauya[j] = dtauy[j] + ctauya[j];

            /* 累积风速三次方（或耦合器所需的相关摩擦诊断量）。 */
            cust3a[j] = dust3[j] + cust3a[j];

            /* 累积感热通量及其关于表面温度的线性化系数。 */
            cshfla[j] = dshfl[j] + cshfla[j];
            cshdta[j] = dshdt[j] + cshdta[j];

            /* 累积潜热通量及其关于表面温度的线性化系数。 */
            clhfla[j] = dlhfl[j] + clhfla[j];
            clhdta[j] = dlhdt[j] + clhdta[j];

            /* 累积表面层短波和长波辐射通量。 */
            cswfla[j] = dswfl[surface] + cswfla[j];
            clwfla[j] = dlwfl[surface] + clwfla[j];
        }
    }

    /* 当前调用已经完成一次大气通量采样，因此累计次数加一。 */
    *naccua += 1;

    /*
     * 判断当前是否真正执行海气/海冰耦合：
     *
     * 任一条件不满足时，只保留累计结果，等待下一次调用继续累积。
     */
    if (nstep % ncpl_atmos_ice != 0 || nkits != 0)
    {
        return 0;
    }

    const float count = (float)(*naccua);

    /*
     * 到达耦合时刻后，将每个海洋格点上的累计量除以采样次数，
     * 从累计和转换成耦合周期内的算术平均值。
     *
     * 这里使用同一个 count，是因为所有累计数组都对应同一组大气时间步；
     * 本函数只负责完成平均，不负责在平均后清零。
     */
    for (int32_t j = 0; j < nhor; ++j)
    {
        if (dls[j] < 0.5f)
        {
            cheata[j] = cheata[j] / count;
            cpmea[j] = cpmea[j] / count;
            cprsa[j] = cprsa[j] / count;
            ctauxa[j] = ctauxa[j] / count;
            ctauya[j] = ctauya[j] / count;
            cust3a[j] = cust3a[j] / count;
            cshfla[j] = cshfla[j] / count;
            cshdta[j] = cshdta[j] / count;
            clhfla[j] = clhfla[j] / count;
            clhdta[j] = clhdta[j] / count;
            cswfla[j] = cswfla[j] / count;
            clwfla[j] = clwfla[j] / count;
        }
    }

    /*
     * 返回 1，通知调用方当前耦合输入已经准备完成。
     */
    return 1;
}

/*
 * 完成一次海气/海冰耦合后的收尾工作。
 * 主要包含两部分：
 *   1. 清零所有大气通量累计数组，为下一个耦合周期做准备；
 *   2. 对海洋格点，将耦合器返回的新海洋和海冰状态复制回大气侧。
 * 最后将累计计数器 naccua 重置为 0。
 */
void sea_finish_coupling(
    int32_t nhor, int32_t *naccua, const float *dls,
    const float *cts, const float *csst, const float *cmld,
    const float *cicec, const float *ciced, const float *csnow,
    const float *csmelt, const float *csndch,
    float *dts, float *dsst, float *dmld, float *dicec, float *diced,
    float *dsnow, float *dsmelt, float *dsndch,
    float *cheata, float *cpmea, float *cprsa, float *ctauxa,
    float *ctauya, float *cust3a, float *cshfla,
    float *cshdta, float *clhfla, float *clhdta, float *cswfla,
    float *clwfla)
{
    /*
     * 遍历全部水平格点。
     * 累计数组对所有格点都清零；
     */
    for (int32_t j = 0; j < nhor; ++j)
    {
        /*
         * 清空上一个耦合周期的大气通量累计值。
         */
        cheata[j] = 0.0f;
        cpmea[j] = 0.0f;
        cprsa[j] = 0.0f;
        ctauxa[j] = 0.0f;
        ctauya[j] = 0.0f;
        cust3a[j] = 0.0f;
        cshfla[j] = 0.0f;
        cshdta[j] = 0.0f;
        clhfla[j] = 0.0f;
        clhdta[j] = 0.0f;
        cswfla[j] = 0.0f;
        clwfla[j] = 0.0f;

        /*
         * 仅对海洋格点读取本次耦合后更新的海洋和海冰状态。
         */
        if (dls[j] < 0.5f)
        {
            dsst[j] = csst[j];   /* 更新海表温度 */
            dts[j] = cts[j];     /* 更新表面温度 */
            dmld[j] = cmld[j];   /* 更新混合层深度 */
            dicec[j] = cicec[j]; /* 更新海冰覆盖率或浓度 */
            diced[j] = ciced[j]; /* 更新海冰厚度 */
            dsnow[j] = csnow[j]; /* 更新海冰积雪状态 */
            dsmelt[j] = csmelt[j];
            dsndch[j] = csndch[j];
        }
    }

    /*
     * 当前耦合周期已经结束，累计次数重置为零。
     */
    *naccua = 0;
}

/*
 * 根据最新海洋/海冰状态，更新大气模式使用的表面边界条件。
 *
 * 仅对 dls[j] < 0.5 的海洋格点执行更新，包括：
 *
 *   1. 将海表温度写入大气温度场表面层；
 *   2. 根据表面温度和局地表面气压计算饱和比湿；
 *   3. 更新大气表面层比湿；
 *   4. 根据海冰覆盖率计算湿度参数；
 *   5. 根据海冰覆盖率计算反照率；
 *   6. 根据 Charnock 关系计算开阔海面粗糙度；
 *   7. 将开阔海面粗糙度与海冰粗糙度线性混合。
 */
void sea_update_surface(
    int32_t nhor, int32_t surface_level,
    float rdbrv, float ra1, float ra2, float ra4, float tmelt,
    float charnock, float gascon, float gravity,
    float albsea, float albice, float dz0sea, float dz0ice,
    float dwetfacsea, float dwetfacice,
    const float *dls, const float *dts, const float *dicec,
    float *dt, float *dq, const float *dp,
    const float *dtaux, const float *dtauy,
    float *dqs, float *dwetfac, float *dalb, float *dz0)
{
    /*
     * 湿度转换修正系数，与初始化函数中的定义一致。
     */
    const float humidity_correction = 1.0f / rdbrv - 1.0f;

    for (int32_t j = 0; j < nhor; ++j)
    {
        /* 仅更新海洋格点的表面边界条件。 */
        if (dls[j] < 0.5f)
        {
            /*
             * 当前水平格点在大气温度/湿度三维场表面层中的索引。
             */
            const size_t surface = field_index(j, surface_level, nhor);

            /*
             * 将海洋表面温度写入大气温度场的表面层。
             */
            dt[surface] = dts[j];

            const float saturation_mixing_ratio =
                rdbrv * ra1 *
                expf(ra2 * (dt[surface] - tmelt) /
                     (dt[surface] - ra4)) /
                dp[j];
            const float saturation_specific_humidity =
                saturation_mixing_ratio /
                (1.0f - humidity_correction * saturation_mixing_ratio);

            /*
             * 更新二维表面饱和比湿和三维大气表面层比湿。
             */
            dqs[j] = saturation_specific_humidity;
            dq[surface] = saturation_specific_humidity;

            /*
             * 计算表面应力矢量的模：
             *     stress = sqrt(taux^2 + tauy^2)
             */
            const float stress = sqrtf(dtaux[j] * dtaux[j] + dtauy[j] * dtauy[j]);

            /*
             * 根据 Charnock 关系计算动态海面粗糙度。
             */
            const float charnock_roughness =
                charnock * stress * gascon * dt[surface] /
                (gravity * dp[j]);

            /*
             * 对动态计算的粗糙度设置最低值。
             * 即使应力很小，开阔海面的粗糙度也不会低于 dz0sea。
             */
            const float open_ocean_roughness =
                fmaxf(charnock_roughness, dz0sea);

            /*
             * 根据海冰覆盖率，对海面和海冰的湿度参数做线性混合。
             */
            dwetfac[j] = dwetfacsea * (1.0f - dicec[j]) + dwetfacice * dicec[j];

            /*
             * 根据海冰覆盖率计算格点平均反照率。
             * 海冰反照率随温度降低而增加，但最大不超过 albice。
             */
            dalb[j] = albsea * (1.0f - dicec[j]) + dicec[j] * fminf(albice, 0.5f + 0.025f * (273.0f - dts[j]));

            /*
             * 根据海冰覆盖率，对动态海面粗糙度和固定海冰粗糙度
             * 进行线性混合：
             *   dicec = 0：使用开阔海面粗糙度 open_ocean_roughness
             *   dicec = 1：使用海冰粗糙度 dz0ice
             */
            dz0[j] = open_ocean_roughness * (1.0f - dicec[j]) +
                     dz0ice * dicec[j];
        }
    }
}
