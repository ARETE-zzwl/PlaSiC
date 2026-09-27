#include "linear_algebra.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/*
 * LAPACK 原生采用列主序（column-major）存储矩阵。
 * 调用者传入的 matrix / inverse_matrix 必须与这里期望的
 * LAPACK 内存布局保持一致，否则得到的结果可能相当于对
 * 转置矩阵进行了计算。
 */

/*
 * SGETRF:
 * 对单精度实矩阵执行 LU 分解，并进行部分主元选取。
 *
 * 数学上得到：
 *
 *     A = P * L * U
 *
 * 更准确地说，LAPACK 使用 permutation/pivot 数组描述行交换。
 *
 * 参数：
 *   row_count         - 矩阵行数 M
 *   column_count      - 矩阵列数 N
 *   matrix            - 输入矩阵，同时也是输出缓冲区；
 *                       调用完成后不再保存原始矩阵，而是保存
 *                       L 和 U 的紧凑表示
 *   leading_dimension - 矩阵的 leading dimension
 *   permutation       - pivot 数组，记录 LU 分解过程中发生的行交换
 *   info              - LAPACK 返回状态
 *
 * info:
 *   = 0  : 成功
 *   < 0  : 第 -info 个参数非法
 *   > 0  : U(info, info) 为 0，矩阵奇异，无法正常求逆
 */
extern void sgetrf_(
    const int32_t *row_count, const int32_t *column_count, float *matrix,
    const int32_t *leading_dimension, int32_t *permutation, int32_t *info);

/*
 * SGETRI:
 * 根据 SGETRF 得到的 LU 分解结果计算矩阵逆。
 *
 * 该函数不是直接接收原始矩阵，而是要求：
 *
 *   1. matrix 已经经过 SGETRF 分解；
 *   2. permutation 是对应 SGETRF 产生的 pivot 数组。
 *
 * 参数：
 *   size              - 方阵阶数 N
 *   matrix            - 输入时保存 SGETRF 的 LU 分解结果；
 *                       输出时被覆盖为原矩阵的逆矩阵
 *   leading_dimension
 *   permutation       - SGETRF 返回的 pivot 数组
 *   work              - LAPACK 工作区
 *   work_size         - 工作区长度 LWORK
 *   info              - LAPACK 返回状态
 *
 * info:
 *   = 0  : 成功
 *   < 0  : 第 -info 个参数非法
 *   > 0  : U(info, info) 为 0，矩阵奇异，无法求逆
 */
extern void sgetri_(
    const int32_t *size, float *matrix, const int32_t *leading_dimension,
    const int32_t *permutation, float *work, const int32_t *work_size,
    int32_t *info);

/*
 * SPTSV solves A X = B for a real, symmetric positive-definite
 * tridiagonal matrix A.  The diagonal, off-diagonal and right-hand side
 * are overwritten by LAPACK; on success right_hand_side contains X.
 */
extern void sptsv_(
    const int32_t *size, const int32_t *right_hand_side_count,
    float *diagonal, float *off_diagonal, float *right_hand_side,
    const int32_t *leading_dimension, int32_t *info);

/*
 * 检查 LAPACK 调用结果。
 *
 * 当前实现采用 fail-fast 策略：
 * 只要 LAPACK 返回的 info != 0，就直接调用 abort() 终止进程。
 */
static void abort_on_lapack_error(int32_t info)
{
    if (info != 0)
    {
        abort();
    }
}

/*
 * 计算 size x size 单精度浮点方阵的逆矩阵。
 * 算法流程：
 *   1. 为 LAPACK pivot 数组分配 permutation；
 *   2. 为 SGETRI 分配工作区 work；
 *   3. 如果输入输出不是同一缓冲区，将 matrix 拷贝到 inverse_matrix；
 *   4. 调用 SGETRF：
 *          A -> LU
 *   5. 调用 SGETRI：
 *          LU -> A^{-1}
 *   6. 释放临时内存。
 * 注意：
 *   当前接口没有返回错误状态。
 *   LAPACK 出错或内存分配失败时，函数会直接 abort()。
 */
void matrix_inverse(
    int32_t size, const float *matrix, float *inverse_matrix)
{
    /*
     * LAPACK SGETRF 返回的 pivot / permutation 数组。
     *
     * 数组长度至少需要 min(M, N)。
     * 由于这里处理的是 size x size 方阵，因此分配 size 个 int32_t。
     */
    int32_t *permutation;

    /*
     * SGETRI 使用的临时工作区。
     *
     * LAPACK 的 SGETRI 允许通过 workspace query 获取性能更优的
     * LWORK；当前实现直接使用 work_size = size，即使用合法的
     * 最小规模工作区，而没有查询最优工作区大小。
     */
    float *work;

    /*
     * LAPACK 状态码。
     */
    int32_t info = 0;

    /*
     * SGETRI 工作区长度。
     * 当前固定设置为 size。
     * 对 SGETRI 来说：
     *   LWORK >= max(1, N)
     * 即可满足基本要求。
     * 更大的、经过 LAPACK workspace query 获得的工作区
     * 可能具有更好的性能，但不影响算法正确性。
     */
    const int32_t work_size = size;

    /*
     * 矩阵元素总数：
     */
    size_t value_count;

    if (size <= 0)
    {
        return;
    }

    /*
     * 计算矩阵元素总数。
     */
    value_count = (size_t)size * (size_t)size;

    /*
     * 分配 LU 分解需要的 pivot 数组。
     * 对 N x N 方阵，需要 N 个 LAPACK INTEGER。
     */
    permutation = (int32_t *)malloc((size_t)size * sizeof(int32_t));

    /*
     * 分配 SGETRI 工作区。
     * work_size 当前等于 size，因此分配 size 个 float。
     */
    work = (float *)malloc((size_t)work_size * sizeof(float));

    /*
     * 任意一次内存分配失败都视为致命错误。
     */
    if (permutation == NULL || work == NULL)
    {
        free(work);
        free(permutation);
        abort();
    }

    if (matrix != inverse_matrix)
    {
        memcpy(inverse_matrix, matrix, value_count * sizeof(float));
    }

    /*
     * LAPACK 在传入的 matrix buffer 上原地完成分解和求逆。
     * 第一步：LU 分解。
     * 调用前：
     *     inverse_matrix = A
     * 调用后：
     *     inverse_matrix = L/U 的紧凑存储形式
     *     permutation    = 行主元交换信息
     * 参数说明：
     *   &size           -> M = N = size
     *   inverse_matrix  -> 待分解矩阵
     *   &size           -> LDA = size
     *   permutation     -> pivot 数组
     *   &info           -> 返回状态
     */
    sgetrf_(&size, &size, inverse_matrix, &size, permutation, &info);

    abort_on_lapack_error(info);

    /*
     * 第二步：根据 LU 分解计算逆矩阵。
     */
    sgetri_(
        &size, inverse_matrix, &size, permutation, work, &work_size, &info);

    abort_on_lapack_error(info);

    free(work);
    free(permutation);
}

/**
 * 求解对称正定三对角线性方程组 A * x = b
 * 该函数调用 LAPACK 中的 SPTSV 例程进行求解。
 *
 * size              矩阵 A 的阶数，即方程组的未知数个数
 * diagonal          三对角矩阵 A 的主对角线元素，长度为 size。注意：SPTSV 调用后，该数组内容可能会被修改
 * off_diagonal      三对角矩阵 A 的次对角线元素，长度为 size - 1。由于 A 为对称矩阵，上对角线和下对角线相同。注意：SPTSV 调用后，该数组内容可能会被修改
 * right_hand_side   方程组右端项 b，长度为 size。求解完成后，该数组会被原地覆盖为解向量 x
 *
 * @return 0    求解成功
 * @return -1   输入的 size 非法（size <= 0）
 * @return >0   LAPACK SPTSV 返回的错误信息，通常表示矩阵无法完成正定分解
 */
int32_t solve_spd_tridiagonal(
    int32_t size, float *diagonal, float *off_diagonal,
    float *right_hand_side)
{
    // 右端项的个数。
    // 当前函数只求解一个方程组 A * x = b，因此设置为 1。
    const int32_t right_hand_side_count = 1;

    // LAPACK 返回的状态码：
    int32_t info = 0;

    // 矩阵阶数必须大于 0。
    if (size <= 0)
    {
        return -1;
    }

    // 调用 LAPACK 的 SPTSV 求解对称正定三对角方程组：
    //
    //     A * X = B
    sptsv_(&size, &right_hand_side_count, diagonal, off_diagonal,
           right_hand_side, &size, &info);

    // 将 LAPACK 的状态码返回给调用者。
    return info;
}
