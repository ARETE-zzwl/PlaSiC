#include "gaussmod_kernels.h"

#include <math.h>
#include <stddef.h>

/*
 * 计算与 Gauss-Legendre 求积相关的归一化 Legendre 多项式值。
 *
 * 参数：
 *   degree : 多项式阶数 n
 *   point  : 计算点 x，通常位于 [-1, 1]
 *
 * 返回值：
 *   degree 阶归一化 Legendre 多项式在 point 处的值。
 */
double legendre_norm(int32_t degree, double point)
{
    /*
     * 将 x 转换为角度：
     *     x = cos(angle)
     */
    double angle = acos(point);

    /*
     * 初始最高阶项系数设为 1。
     */
    double factor = 1.0;

    /* 保存余弦级数累加结果。 */
    double sum = 0.0;

    /*
     * 保存当前一次循环计算出来的项。
     * 注意：
     *   循环结束后，偶数 degree 的情况下还会再次使用 term，
     *   因此 term 定义在循环外部。
     */
    double term = 0.0;

    /* 用于生成下一项 factor 的中间递推系数。 */
    double coefficient;

    /* 当前展开项对应的阶数。 */
    int32_t index;

    /*
     * 从 degree 开始，每次下降 2：
     *
     *   degree, degree-2, degree-4, ...
     */
    for (index = degree; index >= 0; index -= 2)
    {
        term = factor * cos(angle * (double)index);

        sum = sum + term;

        /*
         * 计算下一项所需的递推系数：
         */
        coefficient =
            (double)((degree - index + 1) * (degree + index)) * 0.5;

        /*
         * 更新下一项的系数 factor。
         */
        factor = factor * coefficient /
                 (coefficient + (double)(index - 1));
    }

    /*
     * 如果 degree 为偶数，最后一项对应 index == 0：
     *
     *     cos(0 * angle) = 1
     *
     * 常数项在这种余弦展开中需要取一半权重，因此这里减掉
     * 0.5 * term。
     *
     * 因为循环中已经加入了一整个 term，所以执行：
     *
     *     sum -= 0.5 * term
     *
     * 等价于最终只保留 0.5 * term。
     */
    if (degree % 2 == 0)
    {
        sum = sum - 0.5 * term;
    }

    /*
     * 接下来计算 degree 相关的归一化因子。
     *
     * angle 变量在这里被复用。
     *
     * 初始值为 sqrt(2)。
     */
    angle = sqrt(2.0);

    for (index = 1; index <= degree; ++index)
    {
        angle = angle *
                sqrt(1.0 - 0.25 / (double)(index * index));
    }

    /*
     * 返回：
     *     归一化系数 * 余弦级数结果
     */
    return angle * sum;
}

/*
 * 计算归一化 Legendre 多项式导数的倒数。
 *
 * 参数：
 *   degree : Legendre 多项式阶数 n
 *   point  : 当前迭代点 x
 *
 * 返回值：
 *   1 / Q_n'(x)，其中 Q_n 为 legendre_norm() 给出的
 *   归一化 Legendre 多项式。
 *
 *   gauss_inigau() 中与 legendre_norm() 配合使用：
 *
 *       z3 = legendre_norm(n, x) *
 *            legendre_norm_deriv_reciprocal(n, x);
 *       x  = x - z3;
 *
 *   二者之积构成 Newton 法的根修正步长 Q_n(x) / Q_n'(x)。
 *
 * 说明：
 *   本函数仅供本文件内部使用，因此声明为 static。
 */
static double legendre_norm_deriv_reciprocal(
    int32_t degree, double point)
{
    /*
     * 构造由 degree 阶和 degree-1 阶归一化 Legendre
     * 多项式组成的中间量 z。
     */
    const double z =
        point * legendre_norm(degree, point) -
        sqrt(((double)(degree + degree) + 1.0) /
             ((double)(degree + degree) - 1.0)) *
            legendre_norm(degree - 1, point);

    /*
     * 返回 Newton 修正中使用的比例因子：
     *
     *              x^2 - 1
     *     ---------------------------
     *              degree * z
     */
    return (point * point - 1.0) / ((double)degree * z);
}

/*
 * 初始化 Gaussian / Gauss-Legendre 求积所需要的：
 *
 *   1. abscissas : 求积节点（Legendre 多项式的零点）
 *   2. weights   : 对应节点的 Gaussian quadrature 权重
 *
 * 参数：
 *   latitudes  : 节点总数，同时也是 Legendre 多项式的阶数
 *   abscissas  : 输出数组，长度至少为 latitudes
 *   weights    : 输出数组，长度至少为 latitudes
 *
 * 算法概要：
 *
 *   1. 利用 Legendre 多项式零点的渐近公式获得初值；
 *   2. 使用 Newton 迭代求出正半轴上的零点；
 *   3. 利用 Legendre 多项式零点关于 0 对称的性质，
 *      同时生成对应负半轴零点；
 *   4. 计算每个节点对应的 Gaussian quadrature 权重。
 *
 * 对称性：
 *     x_i = -x_{N-1-i}
 *     w_i =  w_{N-1-i}
 */
void gauss_inigau(
    int32_t latitudes, double *abscissas, double *weights)
{
    /*
     * Newton 方法允许的最大迭代次数。
     * 一般 Legendre 根的初始估计已经很好，通常远不到
     * 50 次就能够收敛。
     */
    static const int32_t iterations = 50;

    /*
     * π 的硬编码近似值。
     * 用于构造 Legendre 零点的初始角度估计。
     */
    static const double pi = 3.14159265358979;

    /*
     * Newton 迭代停止阈值。
     * 当一次 Newton 修正量 |z3| < 1e-16 时，
     * 认为零点已经收敛。
     * 该阈值已经接近 double 精度的机器精度数量级。
     */
    static const double epsilon = 1.0e-16;

    /*
     * 以下 z0 ~ z5 均为算法中的临时变量。
     */
    double z0;
    double z1;
    double z2;
    double z3;
    double z4;
    double z5;

    int32_t latitude;

    int32_t iteration;

    if (latitudes <= 0 || abscissas == NULL || weights == NULL)
    {
        return;
    }

    /*
     * 零点初值公式中使用的公共角度：
     *     z0 = π / (2N + 1)
     */
    z0 = pi / (double)(2 * latitudes + 1);

    /*
     * 初值的高阶修正参数：
     *     z1 = 1 / (8N²)
     */
    z1 = 1.0 / (double)(latitudes * latitudes * 8);

    /*
     * Gaussian 权重公式中的公共因子：
     *     z4 = 2 / N²
     */
    z4 = 2.0 / (double)(latitudes * latitudes);

    /*
     * Legendre 多项式的根关于 0 对称，因此这里只需求
     * 正半轴上的前 N/2 个根。
     * 每求出一个正根 z2，就同时写入：
     *     +z2
     *     -z2
     * 两个对称节点。
     */
    for (latitude = 1; latitude <= latitudes / 2; ++latitude)
    {
        /*
         * 构造第 latitude 个 Legendre 根的初始角度：
         *     theta ≈ π/(2N+1) * (2*latitude - 1/2)
         * 这是基于 Legendre 零点渐近分布得到的初值。
         */
        z2 = z0 * ((double)(2 * latitude) - 0.5);

        /*
         * 将角度估计转换成 x = cos(theta)。
         * 同时增加：
         *     z1 / tan(z2)
         * 对初始角度进行一个有限阶修正，从而让 Newton
         * 方法的初值更加接近真实零点。
         */
        z2 = cos(z2 + z1 / tan(z2));

        /*
         * 使用 Newton 方法迭代求解：
         *     Q_N(z2) = 0
         * 其中 legendre_norm() 给出归一化的 N 阶
         * Legendre 多项式。
         */
        for (iteration = 1; iteration <= iterations; ++iteration)
        {
            /*
             * 计算 Newton 修正量。
             * 从后面的：
             *     z2 = z2 - z3
             * 可以看到 z3 对应 Newton 法中的：
             *     f(x) / f'(x)
             * 这里通过 legendre_norm() 与
             * legendre_norm_deriv_reciprocal() 的组合
             * 来得到该修正量。
             */
            z3 = legendre_norm(latitudes, z2) *
                 legendre_norm_deriv_reciprocal(latitudes, z2);

            /*
             * Newton 更新：
             *     x_{k+1} = x_k - Δx
             */
            z2 = z2 - z3;

            /*
             * 如果本次位置修正已经足够小，则认为收敛。
             */
            if (fabs(z3) < epsilon)
            {
                break;
            }
        }

        /*
         * 在最终求出的 N 阶 Legendre 零点 z2 上，
         * 计算 N-1 阶归一化 Legendre 多项式。
         * 再除以 sqrt(N - 1/2)，形成后续权重公式所需的
         * 辅助量 z5。
         */
        z5 = legendre_norm(latitudes - 1, z2) /
             sqrt((double)latitudes - 0.5);

        /*
         * 保存当前正半轴上的零点。
         * latitude 从 1 开始，而 C 数组下标从 0 开始，
         * 所以写入 latitude - 1。
         */
        abscissas[latitude - 1] = z2;

        /*
         * 计算当前 Gauss-Legendre 节点对应的求积权重：
         *              2        1 - z2²
         *     w = ----------- * --------
         *             N²          z5²
         */
        weights[latitude - 1] =
            z4 * (1.0 - z2 * z2) / (z5 * z5);

        /*
         * 利用 Legendre 多项式零点的对称性：
         *     如果 z2 是一个根，则 -z2 也是一个根。
         * 将对应负半轴节点写入数组另一端。
         */
        abscissas[latitudes - latitude] = -z2;

        /*
         * Gauss-Legendre 求积权重同样关于 0 对称：
         *     w(+x) = w(-x)
         */
        weights[latitudes - latitude] = weights[latitude - 1];
    }
}
