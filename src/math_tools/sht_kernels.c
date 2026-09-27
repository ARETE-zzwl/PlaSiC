/*
 * 本文件是 PlaSiC 与 SHTns 球谐变换库之间的接口层。
 *
 * 主要功能：
 * 1. 管理 SHTns 变换状态；
 * 2. 完成 PlaSiC 数据格式与 SHTns 内部格式之间的转换；
 * 3. 提供标量场、矢量场、散度、涡度以及趋势项计算接口。
 */

#include "sht_kernels.h"
#include "mpi/plasic_mpi.h"

#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <shtns.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * 球谐变换状态结构体。
 *
 * 用于保存一次 SHTns 计算生命周期内需要共享的信息：
 * - config: SHTns 核心配置；
 * - 网格尺寸及 MPI 分区信息；
 * - 谱模式数量；
 * - PlaSiC 与 SHTns 归一化转换系数；
 * - 空间域和谱域临时缓存。
 */
typedef struct sht_transform_state
{
    /*
     * SHTns 核心配置对象。
     * config != NULL 表示当前模块已经完成初始化。
     */
    shtns_cfg config;

    int nlat;
    int nlon;

    /*
     * 当前 MPI rank 在全球纬度方向上的起始偏移。
     *
     * PlaSiC 使用纬向分块并行：
     *
     * 全球纬度:
     *   | rank0 | rank1 | rank2 | ...
     *
     * latitude_offset 表示当前 rank 对应的第一条纬圈编号。
     */
    int latitude_offset;

    /*
     * 当前 MPI rank 持有的局地网格点数量。
     *
     * 通常为：
     *
     * local_points = nlpp(number of latitude points per processor) * nlon
     *
     * 其中 nlpp 为当前进程负责的纬度层数。
     */
    size_t local_points;

    /*
     * 全球网格点总数量。
     * 即：
     * global_points = nlat * nlon
     */
    size_t global_points;

    /*
     * 球谐谱空间中的模式数量。
     * 对应 SHTns 中的 nlm。
     * 每一个 mode 对应一个球谐基函数：
     * Y_l^m
     * 模式按照 SHTns 内部顺序排列。
     */
    size_t spectral_modes;

    /*
     * PlaSiC 与 SHTns 球谐系数归一化之间的转换因子。
     * SHTns:
     *   使用 globally orthonormal spherical harmonics
     * PlaSiC:
     *   使用不同的 Legendre 函数归一化方式
     * 因此需要通过该系数进行谱系数转换。
     */
    double coefficient_scale;

    /*
     * 全球网格临时缓存 A。
     * 用途：
     * - MPI gather 后保存完整全球场；
     * - 标量场转换时作为中间存储。
     * float 类型保持与 PlaSiC 输入输出一致。
     */
    float *global_a;

    /*
     * 全球网格临时缓存 B。
     * 主要用于：
     * - 矢量场计算中的第二个分量；
     * - 保存第二个全球网格变量。
     */
    float *global_b;

    /*
     * 临时保存转换后的 PlaSiC 格式谱系数 A。
     * 存储形式：
     * [Re(mode0), Im(mode0),
     *  Re(mode1), Im(mode1),
     *  ...]
     * 用于趋势项等计算过程中的中间结果。
     */
    float *temporary_spectral_a;

    /*
     * 临时保存转换后的 PlaSiC 格式谱系数 B。
     * 与 temporary_spectral_a 对应，
     * 用于矢量场或多算子组合计算。
     */
    float *temporary_spectral_b;

    /*
     * SHTns 空间域缓存 A。
     * 注意：
     * SHTns 内部要求 double 精度空间数组。
     * 用于：
     * - spat_to_SH()
     * - SH_to_spat()
     * 等标量球谐变换接口。
     */
    double *grid_a;

    /*
     * SHTns 空间域缓存 B。
     * 用于矢量球谐变换：
     * spat_to_SHsphtor()
     * 中的第二个矢量分量。
     */
    double *grid_b;

    /*
     * SHTns 谱空间缓存 A。
     * complex 类型：
     * spectral_a = Re + i * Im
     * 保存一个球谐模式对应的复数系数。
     */
    cplx *spectral_a;

    /*
     * SHTns 谱空间缓存 B。
     * 用于矢量球谐计算中的第二组谱系数。
     */
    cplx *spectral_b;
} sht_transform_state;

static sht_transform_state transform;

static int is_initialized(void)
{
    return transform.config != NULL;
}

static int valid_fields(int levels, const float *first, const float *second)
{
    return is_initialized() && levels >= 0 && first != NULL && second != NULL;
}

static size_t spectral_values(void)
{
    return 2U * transform.spectral_modes;
} 
// 2U 无符号整数（unsigned int）的常量 2

/*
 * 网格数据收集函数。
 *
 * 输入：
 *   local  - 当前进程负责的局地网格。
 *
 * 输出：
 *   global - 拼接后的全球网格。
 *
 * MPI 开启时执行并行 gather；
 * 非 MPI 情况下直接复制数据。
 */
static int gather_grid(const float *local, float *global)
{
    const int32_t status = mp_allgather_grid(
        global, local, (int32_t)transform.global_points,
        (int32_t)transform.local_points, 1);
    /* 1 表示 1 个层次。 */

    return status == 0 ? PLASIC_SHT_OK : PLASIC_SHT_MPI_FAILED;
}

/*
 * 将 PlaSiC 的谱数组转换为 SHTns 复数谱。
 *
 * PlaSiC:
 *   [实部, 虚部, 实部, 虚部, ...]
 *
 * SHTns:
 *   每个球谐模式对应一个 complex 类型变量。
 */
static void unpack_scalar_spectrum(const float *plasic)
{
    size_t mode;

    // 遍历所有谱模态，解包复数谱系数
    for (mode = 0; mode < transform.spectral_modes; ++mode)
    {
        // 零频模态无虚部，其余模态读取虚部
        const double imaginary = transform.config->mi[mode] == 0U
                                     ? 0.0
                                     : (double)plasic[2U * mode + 1U];

        // 从交错存储的实部/虚部数组中恢复复数谱，并应用缩放系数
        transform.spectral_a[mode] =
            transform.coefficient_scale *
            ((double)plasic[2U * mode] + I * imaginary);
    }
}
/*
 * 将 SHTns 复数谱重新打包为 PlaSiC 使用的交错实数数组。
 */
static void pack_scalar_spectrum(float *plasic)
{
    size_t mode;

    // 遍历所有谱模态，将复数谱系数打包为实部/虚部数组
    for (mode = 0; mode < transform.spectral_modes; ++mode)
    {
        // 保存实部，并撤销系数缩放
        plasic[2U * mode] = (float)(creal(transform.spectral_a[mode]) /
                                    transform.coefficient_scale);

        // 有虚部的模态保存虚部，否则零频等特殊模态虚部置零
        plasic[2U * mode + 1U] =
            transform.config->mi[mode] == 0U
                ? 0.0f
                : (float)(cimag(transform.spectral_a[mode]) /
                          transform.coefficient_scale);
    }
}

/*
 * 从一个全局双精度网格数组 global 中提取当前进程/任务负责的局部网格区域，复制到单精度数组 local 中，并根据 sign 调整符号。
 */
static void copy_local_grid(const double *global, float *local, double sign)
{
    // 根据纬度偏移计算全局网格中的起始位置
    const size_t offset =
        (size_t)transform.latitude_offset * (size_t)transform.nlon;

    size_t point;

    // 复制局部网格数据，并进行符号调整和类型转换
    for (point = 0; point < transform.local_points; ++point)
    {
        local[point] = (float)(sign * global[offset + point]);
    }
}


/*
 * 对一个垂直层上的标量格点场执行球谐分析（格点空间 -> 谱空间）。
 *
 * 参数：
 *   local
 *     当前 MPI rank 持有的局地纬度块，按“纬度在外、经度在内”的顺序
 *     连续存储。串行运行时，它就是完整的全球格点场。
 *
 *   spectral
 *     PlaSiC 格式的输出谱数组。每个球谐模态占两个 float，依次存储
 *     实部和虚部；m = 0 模态的虚部由 pack_scalar_spectrum() 置零。
 *
 *   divide_by_cosine_squared
 *     为 0 时，直接分析输入标量 f。
 *
 *     为非零值时，先分析
 *
 *                         f
 *                  ---------------- ,
 *                    cos^2(latitude)
 *
 *     再将结果写入 spectral。这个选项不是 SHTns 高斯求积的权重
 *     修正；SHTns 会自行处理球谐分析所需的求积权重。它只用于撤销
 *     PlaSiC 动量方程中特意预乘的 Robert-form 度量因子。
 *
 *     PlaSiC 保存的水平风为 Robert 形式：
 *
 *         U = u cos(latitude),    V = v cos(latitude),
 *
 *     因而 U^2 + V^2 = cos^2(latitude) (u^2 + v^2)。在 sht_mktend()
 *     中，kinetic_plus_geopotential 被构造成
 *
 *         2 cos^2(latitude) (E + Phi),
 *
 *     其中 E = (u^2 + v^2) / 2。除以 cos^2(latitude) 后得到
 *     2(E + Phi)；随后 sht_mktend() 再乘 0.5*l*(l+1)，从而得到
 *     散度方程所需的 -Laplacian(E + Phi) 的谱系数。
 *
 *     因此，该选项只能用于确实含有上述 cos^2(latitude) 预因子的
 *     特殊输入；普通温度、湿度、地面气压等标量分析必须传 0。
 */
static int analyse_scalar_level(
    const float *local, float *spectral, int divide_by_cosine_squared)
{
    size_t point;

    /*
     * SHTns 在这里对完整的全球网格执行变换。因此 MPI 模式下先把
     * 各 rank 的连续纬度块收集到 transform.global_a；串行模式下
     * gather_grid() 退化为一次直接复制。
     */
    int status = gather_grid(local, transform.global_a);

    if (status != PLASIC_SHT_OK)
    {
        return status;
    }

    /*
     * PlaSiC 的物理场使用 float，SHTns 的空间数组使用 double。
     * 在复制和类型转换的同时完成可选的 Robert-form 度量还原。
     */
    for (point = 0; point < transform.global_points; ++point)
    {
        double value = (double)transform.global_a[point];

        if (divide_by_cosine_squared)
        {
            /*
             * 格点数组按 [latitude][longitude] 展平，所以整除 nlon
             * 即可得到当前格点所属的全球纬度编号。
             */
            const size_t latitude = point / (size_t)transform.nlon;

            /*
             * SHTns 中 st = sin(theta)，theta 是余纬（colatitude）；若 latitude 表示
             * 地理纬度，则 theta = pi/2 - latitude，因此
             *
             *     st = sin(theta) = cos(latitude).
             *
             * Gaussian 纬度不包含两个精确极点，所以 st 不会严格为零。
             */
            const double cosine = transform.config->st[latitude];

            /* 撤销输入场中预先包含的 cos^2(latitude) 度量因子。 */
            value /= cosine * cosine;
        }

        transform.grid_a[point] = value;
    }

    /*
     * SHTns 在内部完成经度 FFT、纬度 Legendre 分析和 Gaussian 求积，
     * 结果按 SHTns 的复数模态排列写入 transform.spectral_a。
     */
    spat_to_SH(transform.config, transform.grid_a, transform.spectral_a);

    /*
     * 转回 PlaSiC 的 [Re(mode0), Im(mode0), Re(mode1), Im(mode1), ...]
     * 布局，同时撤销 PlaSiC 与 SHTns 之间的球谐归一化缩放。
     */
    pack_scalar_spectrum(spectral);

    return PLASIC_SHT_OK;
}


/*
 * 矢量场网格 -> 散度/涡度谱。
 *
 * 输入：
 *   zonal      东西方向分量（指向东为正）；
 *   meridional 南北方向分量（指向北为正）。
 *
 * 输出：
 *   divergence 散度谱；
 *   vorticity  涡度谱。
 *
 * 数学背景：
 *   本接口依赖切向矢量场的 spheroidal/toroidal 势分解：
 *
 *       v = grad_h(S) - r_hat x grad_h(T)
 *
 *   其中 S 为 spheroidal（椭球面）势，决定散度；T 为 toroidal
 *   （环面）势，决定涡度。谱系数与散度/涡度之间仅差 Laplace
 *   特征值 l(l+1)，见本函数末尾的循环。
 *
 * 注意：
 *   调用 spat_to_SHsphtor() 前交换分量顺序并翻转南北符号，
 *   原因是 SHTns 的 theta 轴沿余纬方向（指向南），且参数顺序
 *   固定为 theta 在前——这与 Robert form 无关。SHTns 的
 *   Robert form（u+iv 组合技巧，见 sht_init()）由 SHTns 内部
 *   处理，不影响本层的调用方式。
 */
static int analyse_vector_level(
    const float *zonal, const float *meridional, float *divergence,
    float *vorticity)
{
    size_t point;
    size_t mode;

    // 收集纬向风场（zonal component）到全局网格数组
    int status = gather_grid(zonal, transform.global_a);

    // 检查纬向风场收集是否成功
    if (status != PLASIC_SHT_OK)
    {
        return status;
    }

    // 收集经向风场（meridional component）到全局网格数组
    status = gather_grid(meridional, transform.global_b);

    // 检查经向风场收集是否成功
    if (status != PLASIC_SHT_OK)
    {
        return status;
    }

    /*
     * 坐标分量转换：
     *
     * SHTns 中矢量球谐变换采用球坐标分量：
     *   theta 方向：指向南方（southward）
     *   phi   方向：指向东方（eastward）
     *
     * PlaSiC 中的风场定义与 SHTns 不完全一致：
     *   zonal      -> 纬向风（通常指向东方）
     *   meridional -> 经向风（通常指向北方）
     *
     * 因此需要调整符号和分量顺序，使其符合 SHTns 的输入约定。
     */
    for (point = 0; point < transform.global_points; ++point)
    {
        // SHTns 的 theta 分量指向南方，因此北向 meridional 分量需要取负
        transform.grid_a[point] = -(double)transform.global_b[point];

        // phi 分量与东方向 zonal 分量一致，直接转换
        transform.grid_b[point] = (double)transform.global_a[point];
    }

    /*
     * 矢量球谐分析（正变换）。被调函数的完整全称为：
     *
     *     spat_to_SHsphtor
     *   = spatial -> spherical harmonic, spheroidal-toroidal
     *
     * 即“格点空间切向矢量 -> 球谐 spheroidal/toroidal 谱系数”
     * 的分析变换。输入输出：
     *
     *   输入  grid_a = theta 分量（指向南）= -meridional；
     *         grid_b = phi   分量（指向东）= +zonal；
     *   输出  spectral_a <- S（spheroidal 势，决定散度）；
     *         spectral_b <- T（toroidal   势，决定涡度）。
     *
     * 详细推导见 docs-site 附录
     * 《The Complete Mathematics Behind spat_to_SHsphtor》。
     */
    spat_to_SHsphtor(
        transform.config,
        transform.grid_a,
        transform.grid_b,
        transform.spectral_a,
        transform.spectral_b);

    /*
     * 根据球谐系数计算散度和涡度谱。
     *
     * 对于球谐阶数 l：
     *
     *       Laplace 算子特征值 = l(l+1)
     *
     * 散度谱：
     *       div = -l(l+1) * 对应系数
     *
     * 涡度谱：
     *       vort =  l(l+1) * 对应系数
     *
     * 这里需要除以 coefficient_scale，
     * 因为内部谱系数存储时进行了缩放。
     */
    for (mode = 0; mode < transform.spectral_modes; ++mode)
    {
        // 获取当前球谐模态的阶数 l
        const double degree = (double)transform.config->li[mode];

        // 球面 Laplace 算子的特征值 l(l+1)
        const double eigenvalue = degree * (degree + 1.0);

        // 计算散度谱系数
        const cplx divergence_value =
            -eigenvalue * transform.spectral_a[mode] /
            transform.coefficient_scale;

        // 计算涡度谱系数
        const cplx vorticity_value =
            eigenvalue * transform.spectral_b[mode] /
            transform.coefficient_scale;

        // 保存散度和涡度谱的实部
        divergence[2U * mode] = (float)creal(divergence_value);
        vorticity[2U * mode] = (float)creal(vorticity_value);

        // 判断当前模态是否包含虚部
        if (transform.config->mi[mode] == 0U)
        {
            // 特殊模态（如零阶/对称模态）虚部理论上为零
            divergence[2U * mode + 1U] = 0.0f;
            vorticity[2U * mode + 1U] = 0.0f;
        }
        else
        {
            // 普通复数模态，保存虚部
            divergence[2U * mode + 1U] = (float)cimag(divergence_value);
            vorticity[2U * mode + 1U] = (float)cimag(vorticity_value);
        }
    }

    return PLASIC_SHT_OK;
}

/*
 * 初始化球谐变换模块（整个模块的入口，必须在使用其它接口前调用一次）。
 *
 * 整体流程分为四步：
 * 1. 参数合法性检查；
 * 2. 创建 SHTns 配置对象并绑定网格（含 FFT 计划等预计算）；
 * 3. 把网格/谱的尺寸信息记录到全局状态 transform 中；
 * 4. 分配各临时缓存并做统一的失败清理。
 *
 * 参数说明：
 *
 *   ntru
 *     三角截断（triangular truncation）波数，即最大球谐总波数 l_max，
 *     同时也取 m_max = ntru。例如 T42 对应 ntru = 42。
 *
 *   nlat
 *     全球纬圈（纬度方向）格点数。
 *
 *   nlon
 *     全球经度（经度方向）格点数。
 *
 *   nlpp
 *     number of latitude points per processor：
 *     当前 MPI 进程负责的纬圈数（串行模式下 nlpp = nlat）。
 *
 *   latitude_offset
 *     当前进程在全球纬圈数组中的起始编号（从 0 开始），
 *     用于 MPI gather/scatter 时确定本进程数据在全球场中的位置。
 *
 * 返回值：
 *   PLASIC_SHT_OK                   初始化成功；
 *   PLASIC_SHT_INVALID_ARGUMENT     参数组合不合法；
 *   PLASIC_SHT_INITIALIZATION_FAILED SHTns 创建配置或建网格失败；
 *   PLASIC_SHT_ALLOCATION_FAILED    缓存内存分配失败。
 */
int sht_init(int ntru, int nlat, int nlon, int nlpp, int latitude_offset)
{
    enum shtns_type grid_flags;

    /*
     * 第一步：参数合法性检查。
     *
     * 各条件与物理/数值约束一一对应：
     *   ntru < 1                    截断波数必须至少为 1；
     *   nlpp < 1                    每个进程至少要分到一条纬圈；
     *   latitude_offset < 0         纬度偏移不能为负；
     *   latitude_offset + nlpp > nlat
     *                               本进程负责的纬圈范围不能越出全球纬圈。
     */
    if (ntru < 1 || nlpp < 1 ||
        latitude_offset < 0 || latitude_offset + nlpp > nlat)
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    /*
     * 第二步：创建 SHTns 配置。
     *
     * 先调用 sht_finalize() 清理上一次初始化遗留的资源，
     * 保证本函数可以安全地重复调用。
     */
    sht_finalize();

    /*
     * shtns_create 参数依次为：
     *   lmax = ntru       最大球谐总波数；
     *   mmax = ntri       最大傅里叶阶数（此处同样取 ntru，即三角截断）；
     *   mres = 1          磁场类问题用的 m 共振参数，普通气象场景固定取 1，
     *                     表示 m 取遍 0..mmax 的所有整数；
     *   sht_orthonormal   采用正交归一化（globally orthonormal）的球谐基，
     *                     即 ∫ |Y_l^m|^2 dΩ = 1 的约定。
     */
    transform.config = shtns_create(ntru, ntru, 1, sht_orthonormal);
    if (transform.config == NULL)
    {
        // 配置对象创建失败（通常是内部内存分配出错）
        return PLASIC_SHT_INITIALIZATION_FAILED;
    }

    /*
     * 网格构建选项（按位或组合）：
     *   sht_quick_init     允许 SHTns 使用 FFTW 快速路径和预计算表，
     *                      牺牲少量初始化时间换取变换速度；
     *   SHT_PHI_CONTIGUOUS 内存布局约定：同一纬圈上的经度方向（phi）
     *                      数据在内存中连续，便于 SIMD/FFT 向量化；
     *   SHT_ROBERT_FORM    矢量变换使用 Robert form：
     *                      把风的两个分量组合成 u+iv 复数形式
     *                      参与勒让德变换，一次遍历同时得到
     *                      spheroidal/toroidal 两组谱系数。
     *                      注意：analyse_vector_level() 中的
     *                      南北符号翻转源于 SHTns 的余纬方向
     *                      约定，与该选项无关。
     */
    grid_flags = (enum shtns_type)(sht_quick_init | SHT_PHI_CONTIGUOUS |
                                   SHT_ROBERT_FORM);

    /*
     * 绑定具体网格：
     *   极点偏移 0.0、纬圈数 nlat、经度数 nlon。
     * 这一步会预计算高斯节点、Legendre 多项式表和 FFT 计划，
     * 是初始化中最耗时的一步。
     *
     * 返回值含义：成功时返回谱模式数量 nlm（> 0），
     * 因此 <= 0 一律视为建网格失败；失败时需回收已创建的配置。
     */
    if (shtns_set_grid(transform.config, grid_flags, 0.0, nlat, nlon) <= 0)
    {
        sht_finalize();
        return PLASIC_SHT_INITIALIZATION_FAILED;
    }

    /*
     * 显式开启 Robert form（矢量场以 u+iv / u-iv 组合形式参与变换）。
     * 与 grid_flags 中的 SHT_ROBERT_FORM 呼应，
     * 决定了 analyse_vector_level 里分量符号调整的方式。
     */
    shtns_robert_form(transform.config, 1);
    // robert != 0：启用 Robert form
    // robert == 0：关闭 Robert form

    /*
     * 第三步：把网格与谱的尺寸信息写入全局状态 transform，
     * 后续所有变换接口都从这里读取配置。
     */
    transform.nlat = nlat;
    transform.nlon = nlon;
    transform.latitude_offset = latitude_offset;

    transform.local_points = (size_t)nlpp * (size_t)nlon;

    transform.global_points = (size_t)nlat * (size_t)nlon;

    transform.spectral_modes = (size_t)transform.config->nlm;

    /*
     * 归一化转换因子 sqrt(2π)。
     *
     * SHTns 的正交归一化约定与 PlaSiC 内部使用的 Legendre 函数
     * 归一化相差一个 sqrt(2π) 因子：
     *   PlaSiC -> SHTns：系数乘 coefficient_scale
     *   SHTns -> PlaSiC：系数除 coefficient_scale
     * acos(-1.0) 是避免依赖 M_PI 宏的写法。
     */
    transform.coefficient_scale = sqrt(2.0 * acos(-1.0));

    /*
     * 第四步：分配全部临时缓存。
     *
     * 分两组接口分配：
     * - malloc：PlaSiC 侧缓存，float 单精度，保持与 PlaSiC 数据一致；
     * - shtns_malloc：SHTns 侧缓存，double 双精度且按 SIMD 对齐，
     *   SHTns 变换核函数要求这种对齐方式才能保证性能与正确性。
     */
    transform.global_a =
        (float *)malloc(transform.global_points * sizeof(*transform.global_a));
    transform.global_b =
        (float *)malloc(transform.global_points * sizeof(*transform.global_b));

    transform.temporary_spectral_a =
        (float *)malloc(spectral_values() *
                        sizeof(*transform.temporary_spectral_a));
    transform.temporary_spectral_b =
        (float *)malloc(spectral_values() *
                        sizeof(*transform.temporary_spectral_b));

    transform.grid_a = (double *)shtns_malloc(
        (size_t)transform.config->nspat * sizeof(*transform.grid_a));
    transform.grid_b = (double *)shtns_malloc(
        (size_t)transform.config->nspat * sizeof(*transform.grid_b));

    transform.spectral_a = (cplx *)shtns_malloc(
        transform.spectral_modes * sizeof(*transform.spectral_a));
    transform.spectral_b = (cplx *)shtns_malloc(
        transform.spectral_modes * sizeof(*transform.spectral_b));

    /*
     * 统一检查所有分配是否成功。
     * 只要有一个失败就整体回滚：sht_finalize() 会释放已成功的部分，
     * 并把 transform 清零，使模块回到"未初始化"状态，避免泄漏或半初始化。
     */
    if (transform.global_a == NULL || transform.global_b == NULL ||
        transform.temporary_spectral_a == NULL ||
        transform.temporary_spectral_b == NULL ||
        transform.grid_a == NULL || transform.grid_b == NULL ||
        transform.spectral_a == NULL || transform.spectral_b == NULL)
    {
        sht_finalize();
        return PLASIC_SHT_ALLOCATION_FAILED;
    }
    return PLASIC_SHT_OK;
}

/*
 * 释放球谐模块资源。
 */
void sht_finalize(void)
{
    free(transform.global_a);
    free(transform.global_b);
    free(transform.temporary_spectral_a);
    free(transform.temporary_spectral_b);
    shtns_free(transform.grid_a);
    shtns_free(transform.grid_b);
    shtns_free(transform.spectral_a);
    shtns_free(transform.spectral_b);
    if (transform.config != NULL)
    {
        shtns_destroy(transform.config);
    }
    memset(&transform, 0, sizeof(transform));
    // 将变量 transform 所占用的整块内存全部初始化为 0。
}

/*
 * 标量场：谱空间 -> 网格空间（每层一次逆变换）。
 *
 * 用途：
 *   根据球谐系数恢复物理空间变量，例如把谱空间的地面气压、
 *   温度等场变换回高斯网格，供物理过程参数化或输出使用。
 *
 * 参数：
 *   levels    需要变换的垂直层数；
 *   spectral  输入谱数组，PlaSiC 交错格式：
 *             [Re(mode0), Im(mode0), Re(mode1), Im(mode1), ...]
 *             每层占 spectral_values() 个 float；
 *   grid      输出网格数组，当前进程局地网格，
 *             按“纬度在外、经度在内”连续存储，每层 local_points 个点。
 *
 * 每一层的处理分三步：
 *   1. 解包：PlaSiC 交错实数 -> SHTns 复数谱，并补上归一化缩放；
 *   2. SH_to_spat：SHTns 内部完成合成；
 *   3. 提取本进程纬度块并转回 float 输出。
 */
int sht_scalar_to_grid(int levels, const float *spectral, float *grid)
{
    int level;

    // 入口检查：模块已初始化、层数非负、两个指针均非空
    if (!valid_fields(levels, spectral, grid))
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    for (level = 0; level < levels; ++level)
    {
        /*
         * 第 1 步：解包第 level 层的谱数据。
         * spectral + level*spectral_values() 定位到该层起始位置。
         * 结果写入 transform.spectral_a（SHTns 复数格式）。
         */
        unpack_scalar_spectrum(
            spectral + (size_t)level * spectral_values());

        /*
         * 第 2 步：逆球谐变换（谱 -> 全球网格）。
         * SHTns 直接在全球网格上计算，结果为 double 精度。
         */
        SH_to_spat(transform.config, transform.spectral_a, transform.grid_a);

        /*
         * 第 3 步：从全球网格中截取当前进程负责的纬度块。
         * 符号因子取 +1.0，只做 double->float 转换，不改变数值。
         */
        copy_local_grid(
            transform.grid_a,
            grid + (size_t)level * transform.local_points, 1.0);
    }
    return PLASIC_SHT_OK;
}

int sht_local_scalar_to_grid(
    int levels, const float *local_spectral, int local_spectral_values,
    float *grid)
{
    const size_t global_values = spectral_values();
    /* 局部谱块按 ceil(global / size) 预留；最后一个 rank 可能只填满一部分。 */
    /* Local blocks are reserved as ceil(global / size); the last rank may
     * fill only part of its block. */
    const size_t local_values =
        (global_values + (size_t)mp_size() - 1U) / (size_t)mp_size();
    size_t allocation_values;
    float *global_spectral;
    int status;

    if (!valid_fields(levels, local_spectral, grid) ||
        local_spectral_values < 0 ||
        (size_t)local_spectral_values != local_values ||
        global_values > (size_t)INT32_MAX)
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }
    if (mp_size() == 1)
    {
        return sht_scalar_to_grid(levels, local_spectral, grid);
    }

    if (levels != 0 && global_values > SIZE_MAX / (size_t)levels)
    {
        return PLASIC_SHT_ALLOCATION_FAILED;
    }
    allocation_values = global_values * (size_t)levels;
    if (allocation_values > SIZE_MAX / sizeof(*global_spectral))
    {
        return PLASIC_SHT_ALLOCATION_FAILED;
    }
    global_spectral = (float *)malloc(
        (allocation_values == 0U ? 1U : allocation_values) *
        sizeof(*global_spectral));
    if (global_spectral == NULL)
    {
        return PLASIC_SHT_ALLOCATION_FAILED;
    }

    status = (int)mp_allgather_spectral(
        global_spectral, local_spectral, (int32_t)global_values,
        (int32_t)local_spectral_values, (int32_t)levels);
    if (status == 0)
    {
        status = sht_scalar_to_grid(levels, global_spectral, grid);
    }
    else
    {
        status = PLASIC_SHT_MPI_FAILED;
    }
    free(global_spectral);
    return status;
}

/*
 * 标量场：网格空间 -> 谱空间（每层一次正变换）。
 *
 * 用途：
 *   把高斯网格上的标量场（如温度、比湿、地形高度）分析成
 *   球谐谱系数，是 sht_scalar_to_grid 的逆操作。
 *
 * 参数：
 *   levels    需要变换的垂直层数；
 *   grid      输入网格数组，当前进程局地网格；
 *   spectral  输出谱数组，PlaSiC 交错格式，每层 spectral_values()
 *             个 float。
 *
 * 具体的 gather、精度转换、spat_to_SH 变换和重新打包都在
 * analyse_scalar_level() 中完成；
 * 这里传 0 表示不做除以 cos^2(latitude) 的度量还原（普通标量场不需要）。
 */
int sht_grid_to_scalar(int levels, const float *grid, float *spectral)
{
    int level;

    // 入口检查：模块已初始化、层数非负、两个指针均非空
    if (!valid_fields(levels, grid, spectral))
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    for (level = 0; level < levels; ++level)
    {
        // 对第 level 层执行完整的“格点 -> 谱”分析流程
        const int status = analyse_scalar_level(
            grid + (size_t)level * transform.local_points,
            spectral + (size_t)level * spectral_values(), 0);
        
        if (status != PLASIC_SHT_OK)
        {
            return status;
        }
    }
    return PLASIC_SHT_OK;
}

/*
 * 它把一个标量物理场（例如气压、温度）的球谐谱数据，转换成经纬网格上的水平梯度，并输出北向梯度和东向梯度两个分量。
 *
 * 数学背景：
 *   在球谐空间中，梯度算子的作用非常简单——对每个模态
 *
 *       grad f = (1/a) * [ e_theta * dF/dtheta + e_phi * (im/(sin theta)) F ]
 *
 *   即只涉及对缔合 Legendre 函数的求导和乘 m/(sin theta)，
 *   SHTns 的 SHsph_to_spat() 在逆变换的同时自动完成这些运算。
 *
 * 用途：
 *   计算气压梯度力等需要水平导数的项，输出南北、东西两个分量。
 *
 * 参数：
 *   spectral   输入标量谱，PlaSiC 交错格式；
 *   northward  输出向北分量网格场；
 *   eastward   输出向东分量网格场。
 *
 * 注意：
 *   本接口没有 levels 参数，一次只处理一个水平层
 */
int sht_scalar_to_grid_gradient(
    const float *spectral, float *northward, float *eastward)
{
    // 入口检查：模块已初始化且三个指针均非空
    if (!is_initialized() || spectral == NULL || northward == NULL ||
        eastward == NULL)
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    // 解包 PlaSiC 谱 -> SHTns 复数谱（含归一化缩放）
    unpack_scalar_spectrum(spectral);

    /*
     * 梯度型逆矢量变换：
     * 输入为标量谱 spectral_a，
     * 输出 grid_a = theta 方向（指向南）分量，
     *      grid_b = phi   方向（指向东）分量。
     */
    SHsph_to_spat(
        transform.config, transform.spectral_a,
        transform.grid_a, transform.grid_b);

    /*
     * 坐标方向换算：
     * SHTns 的 theta 分量指向南，而 PlaSiC 需要向北分量，
     * 因此取负号；phi 分量本来就指向东，符号不变。
     */
    copy_local_grid(transform.grid_a, northward, -1.0);
    copy_local_grid(transform.grid_b, eastward, 1.0);
    return PLASIC_SHT_OK;
}

/*
 * 散度/涡度谱 -> 风场网格值（Helmholtz 分解的逆过程）。
 *
 * 参数：
 *   levels            需要处理的垂直层数；
 *   planetary_vorticity
 *                     行星涡度在 (m=0, l=1) 模态上的谱系数大小。
 *                     PlaSiC 谱中存储的是绝对涡度，其行星旋转部分
 *                     恰好落在 (m=0, l=1) 实部上；重建风场前必须先
 *                     扣除它，得到相对涡度后才能反演（见冷启动处
 *                     runtime/runtime_dynamics.c 中的说明）。
 *   divergence        输入散度谱，PlaSiC 交错格式，每层一组；
 *   vorticity         输入涡度谱，PlaSiC 交错格式，每层一组；
 *   zonal_wind        输出纬向风 u 网格场；
 *   meridional_wind   输出经向风 v 网格场（指向北为正）。
 */
int sht_vortdiv_to_wind(
    int levels, float planetary_vorticity, const float *divergence,
    const float *vorticity, float *zonal_wind, float *meridional_wind)
{
    int level;

    // 入口检查：模块已初始化、层数非负、全部指针非空
    if (!is_initialized() || levels < 0 || divergence == NULL ||
        vorticity == NULL || zonal_wind == NULL || meridional_wind == NULL)
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    for (level = 0; level < levels; ++level)
    {
        // 定位第 level 层的散度谱 / 涡度谱起始位置
        const float *level_divergence =
            divergence + (size_t)level * spectral_values();
        const float *level_vorticity =
            vorticity + (size_t)level * spectral_values();
        size_t mode;

        /*
         * 第 1 步：逐模态反解速度势/流函数系数。
         *
         * 由 div = -l(l+1)*S 和 vort = +l(l+1)*T 得：
         *   spectral_a = -div / l(l+1)
         *   spectral_b = +vort / l(l+1)
         * 同时乘 coefficient_scale 完成 PlaSiC -> SHTns 归一化转换。
         */
        for (mode = 0; mode < transform.spectral_modes; ++mode)
        {
            // 当前模态的总波数 l 及 Laplace 特征值 l(l+1)
            const double degree = (double)transform.config->li[mode];
            const double eigenvalue = degree * (degree + 1.0);
            if (eigenvalue == 0.0)
            {
                /*
                 * l = 0 的带状平均模态：除以零无意义。
                 * 常数型的散度/涡度对风场没有贡献，直接置零。
                 */
                transform.spectral_a[mode] = 0.0;
                transform.spectral_b[mode] = 0.0;
            }
            else
            {
                /*
                 * m = 0 模态虚部理论上恒为零，读取时强制置 0，
                 * 避免输入数组中的无效数据污染结果。
                 */
                const double divergence_imaginary =
                    transform.config->mi[mode] == 0U
                        ? 0.0
                        : (double)level_divergence[2U * mode + 1U];
                const double vorticity_imaginary =
                    transform.config->mi[mode] == 0U
                        ? 0.0
                        : (double)level_vorticity[2U * mode + 1U];

                // 涡度实部先复制出来，便于做行星涡度修正
                double vorticity_real =
                    (double)level_vorticity[2U * mode];

                /*
                 * 绝对涡度 -> 相对涡度修正：
                 * 行星旋转（刚体旋转）的涡度只出现在 (m=0, l=1)
                 * 模态实部中，这里减去该分量，使后续反演得到的
                 * 是相对运动的风场。
                 */
                if (transform.config->mi[mode] == 0U &&
                    transform.config->li[mode] == 1U)
                {
                    vorticity_real -= (double)planetary_vorticity;
                }

                // S = -div/(l(l+1))，含归一化缩放
                transform.spectral_a[mode] =
                    -transform.coefficient_scale *
                    ((double)level_divergence[2U * mode] +
                     I * divergence_imaginary) /
                    eigenvalue;

                // T = +vort/(l(l+1))，含归一化缩放
                transform.spectral_b[mode] =
                    transform.coefficient_scale *
                    (vorticity_real + I * vorticity_imaginary) /
                    eigenvalue;
            }
        }

        /*
         * 第 2 步：逆矢量变换（谱 -> 全球网格）。
         * 输出 grid_a = theta 方向（向南）分量，
         *      grid_b = phi   方向（向东）分量。
         */
        SHsphtor_to_spat(
            transform.config, transform.spectral_a, transform.spectral_b,
            transform.grid_a, transform.grid_b);

        /*
         * 第 3 步：方向换算并截取本进程纬度块。
         * 纬向风 u 来自 phi（向东）分量，符号不变；
         * 经向风 v 来自 theta（向南）分量，取负翻转为向北为正。
         */
        copy_local_grid(
            transform.grid_b,
            zonal_wind + (size_t)level * transform.local_points, 1.0);
        copy_local_grid(
            transform.grid_a,
            meridional_wind + (size_t)level * transform.local_points, -1.0);
    }
    return PLASIC_SHT_OK;
}

/*
 * 风场网格值 -> 散度/涡度谱（每层一次矢量正变换）。
 *
 * 用途：
 *   sht_vortdiv_to_wind 的逆操作。把高斯网格上的风场（u 纬向、v 经向）分析成散度谱与涡度谱，供动力学核心计算或诊断输出使用。
 *
 * 参数：
 *   levels            需要处理的垂直层数；
 *   zonal_wind        输入纬向风 u 网格场；
 *   meridional_wind   输入经向风 v 网格场（指向北为正）；
 *   divergence        输出散度谱，PlaSiC 交错格式，每层一组；
 *   vorticity         输出涡度谱，PlaSiC 交错格式，每层一组。
 *
 * 具体的 MPI gather、Robert form 坐标换算、spat_to_SHsphtor
 * 变换以及 l(l+1) 算子作用都在 analyse_vector_level() 中完成。
 */
int sht_wind_to_vortdiv(
    int levels, const float *zonal_wind, const float *meridional_wind,
    float *divergence, float *vorticity)
{
    int level;

    // 入口检查：模块已初始化、层数非负、全部指针非空
    if (!is_initialized() || levels < 0 || zonal_wind == NULL ||
        meridional_wind == NULL || divergence == NULL || vorticity == NULL)
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    for (level = 0; level < levels; ++level)
    {
        // 对第 level 层执行完整的“风场 -> 散度/涡度谱”分析流程
        const int status = analyse_vector_level(
            zonal_wind + (size_t)level * transform.local_points,
            meridional_wind + (size_t)level * transform.local_points,
            divergence + (size_t)level * spectral_values(),
            vorticity + (size_t)level * spectral_values());

        if (status != PLASIC_SHT_OK)
        {
            return status;
        }
    }
    return PLASIC_SHT_OK;
}

/*
 * 标量趋势项（平流 + 源汇）：tendency = 源项 - div(通量)。
 *
 * 数学背景：
 *   标量守恒方程的通量形式：
 *
 *       dq/dt = S - div(q * v)
 *
 *   其中 S 为源项（物理过程贡献），q*v 为标量被风输送的通量。
 *   在谱空间中，两项都是普通的球谐系数数组，可以直接逐元素相减。
 *                     scalar_source
 *                         S
 *                         |
 *                         |
 *                   球谐谱分析
 *                         |
 *                         v
 *                      S_lm
 *                         |
 *                         |
 *                         |         zonal_flux
 *                         |         meridional_flux
 *                         |              |
 *                         |              |
 *                         |       矢量球谐谱分析
 *                         |              |
 *                         |              v
 *                         |         div(F)_lm
 *                         |              |
 *                         |              |
 *                         +------- 减 ----+
 *                                 |
 *                                 v
 *                     tendency_lm
 *                   = S_lm - div(F)_lm
 * 参数：
 *   levels            需要处理的垂直层数；
 *   scalar_source     输入源项网格场 S；
 *   zonal_flux        输入纬向通量分量 (q*u)；
 *   meridional_flux   输入经向通量分量 (q*v)，指向北为正；
 *   tendency          输出趋势谱，PlaSiC 交错格式，每层一组。
 */
int sht_qtend(
    int levels, const float *scalar_source, const float *zonal_flux,
    const float *meridional_flux, float *tendency)
{
    int level;

    // 入口检查：模块已初始化、层数非负、全部指针非空
    if (!is_initialized() || levels < 0 || scalar_source == NULL ||
        zonal_flux == NULL || meridional_flux == NULL || tendency == NULL)
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    for (level = 0; level < levels; ++level)
    {
        // 定位第 level 层的输出谱位置
        float *level_tendency =
            tendency + (size_t)level * spectral_values();
        size_t value;

        /*
         * 第 1 步：分析通量场的散度 div(q*v)。
         *
         * 结果直接写入 level_tendency（复用输出缓冲区，
         * 避免额外的临时分配）；涡度输出写入 temporary_spectral_b，
         * 本接口不需要它，后续会被覆盖。
         */
        int status = analyse_vector_level(
            zonal_flux + (size_t)level * transform.local_points,
            meridional_flux + (size_t)level * transform.local_points,
            level_tendency, transform.temporary_spectral_b);

        if (status != PLASIC_SHT_OK)
        {
            return status;
        }

        /*
         * 第 2 步：分析源项 S 的谱系数。
         * 写入临时缓存 temporary_spectral_a；传 0 表示不做
         * cos^2 度量还原（普通标量源项不需要）。
         */
        status = analyse_scalar_level(
            scalar_source + (size_t)level * transform.local_points,
            transform.temporary_spectral_a, 0);
        if (status != PLASIC_SHT_OK)
        {
            return status;
        }

        /*
         * 第 3 步：合成趋势谱。
         * 注意实部/虚部是交错存储的，因此直接按元素遍历
         * spectral_values() 个 float 即可完成 source - div 运算。
         */
        for (value = 0; value < spectral_values(); ++value)
        {
            level_tendency[value] =
                transform.temporary_spectral_a[value] -
                level_tendency[value];
        }
    }
    return PLASIC_SHT_OK;
}


//                     绝对涡度通量
//                     Fη = (Fx, Fy)
//                          │
//                          │ vector SHT
//                          ▼
//               ┌─────────────────────┐
//               │                     │
//               ▼                     ▼
//        div(Fη) 的谱            curl(Fη) 的谱
//               │                     │
//               │                     └──→ 涡度趋势
//               │
//               │
//               ▼
//            散度趋势
//               ▲
//               │
//      -∇²(E + Φ) 的谱
//               ▲
//               │
//        scalar SHT
//               │
//        2(E + Φ)
//               ▲
//               │
//    2 cos²φ (E + Φ)


// 温度部分：

//        T v = (Tu, Tv)
//               │
//               │ vector SHT
//               ▼
//         div(Tv) 的谱
//               │
//               │ 取负
//               ▼
//               +
//               ▲
//               │
//           S_T 的谱
//               ▲
//               │
//          scalar SHT
//               │
//               ▼

//         ∂T/∂t = S_T - div(Tv)

/*
 * 输入约定（与 dynamical_core/plasic_kernels.c 的构造方式对应）：
 *   kinetic_plus_geopotential 中存放的是
 *       2 * cos^2(latitude) * (E + Phi)，
 *   即预先乘了 Robert form 度量因子 cos^2(latitude)；因此分析时
 *   要传 divide_by_cosine_squared = 1 把它还原成 2(E + Phi)。
 */
int sht_mktend(
    int levels, const float *temperature_source,
    const float *zonal_vorticity_flux,
    const float *meridional_vorticity_flux,
    const float *kinetic_plus_geopotential,
    const float *zonal_temperature_flux,
    const float *meridional_temperature_flux, float *divergence_tendency,
    float *temperature_tendency, float *vorticity_tendency)
{
    int level;

    // 入口检查：模块已初始化、层数非负、全部指针非空
    if (!is_initialized() || levels < 0 || temperature_source == NULL ||
        zonal_vorticity_flux == NULL || meridional_vorticity_flux == NULL ||
        kinetic_plus_geopotential == NULL ||
        zonal_temperature_flux == NULL ||
        meridional_temperature_flux == NULL ||
        divergence_tendency == NULL || temperature_tendency == NULL ||
        vorticity_tendency == NULL)
    {
        return PLASIC_SHT_INVALID_ARGUMENT;
    }

    for (level = 0; level < levels; ++level)
    {
        // 第 level 层在网格数组 / 谱数组中的起始偏移
        const size_t grid_offset =
            (size_t)level * transform.local_points;
        const size_t spectral_offset =
            (size_t)level * spectral_values();

        // 三个输出谱趋势的本层起始位置（后续反复读写）
        float *level_divergence = divergence_tendency + spectral_offset;
        float *level_temperature = temperature_tendency + spectral_offset;
        float *level_vorticity = vorticity_tendency + spectral_offset;
        size_t mode;
        size_t value;

        /*
         * 第 1 步：分析绝对涡度通量矢量。
         *
         * 一次矢量正变换同时得到两个结果：
         *   level_divergence <- div(涡度通量)：散度方程的对流项；
         *   level_vorticity  <- vort(涡度通量)：涡度方程的趋势主体。
         */
        int status = analyse_vector_level(
            zonal_vorticity_flux + grid_offset,
            meridional_vorticity_flux + grid_offset,
            level_divergence, level_vorticity);

        if (status != PLASIC_SHT_OK)
        {
            return status;
        }

        /*
         * 第 2 步：分析动能 + 位势项。
         *
         * 输入含 cos^2 预因子，传 1 触发除以 cos^2(latitude)，
         * 结果 2(E+Phi) 的谱系数暂存于 level_temperature 缓冲区。
         */
        status = analyse_scalar_level(
            kinetic_plus_geopotential + grid_offset,
            level_temperature, 1);
        if (status != PLASIC_SHT_OK)
        {
            return status;
        }

        /*
         * 第 3 步：把能量项贡献累加进散度趋势。
         *
         * 对每个模态乘 0.5*l(l+1) 后加到散度趋势上：
         *   0.5*l(l+1) × [2(E+Phi)] 的谱 = -Laplacian(E+Phi)，
         * 正是散度方程右端所需的气压梯度/能量项。
         * 注意实部、虚部要分别累加。
         */
        for (mode = 0; mode < transform.spectral_modes; ++mode)
        {
            const double degree = (double)transform.config->li[mode];
            const float factor =
                (float)(0.5 * degree * (degree + 1.0));
            level_divergence[2U * mode] +=
                factor * level_temperature[2U * mode];
            level_divergence[2U * mode + 1U] +=
                factor * level_temperature[2U * mode + 1U];
        }

        /*
         * 第 4 步：分析温度平流通量。
         *
         * 与 sht_qtend 相同的模式：
         *   level_temperature <- div(T*v)（复用缓冲区）；
         *   temporary_spectral_b <- 涡度部分，此处不需要。
         */
        status = analyse_vector_level(
            zonal_temperature_flux + grid_offset,
            meridional_temperature_flux + grid_offset,
            level_temperature, transform.temporary_spectral_b);
        if (status != PLASIC_SHT_OK)
        {
            return status;
        }

        // 第 5 步：分析温度源项 S 的谱系数
        status = analyse_scalar_level(
            temperature_source + grid_offset,
            transform.temporary_spectral_a, 0);
        if (status != PLASIC_SHT_OK)
        {
            return status;
        }

        // 第 6 步：温度趋势 = 源项谱 - div(T*v) 谱（逐元素相减）
        for (value = 0; value < spectral_values(); ++value)
        {
            level_temperature[value] =
                transform.temporary_spectral_a[value] -
                level_temperature[value];
        }
    }
    return PLASIC_SHT_OK;
}
