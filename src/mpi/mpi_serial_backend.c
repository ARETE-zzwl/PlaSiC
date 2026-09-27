/*
 * 单进程通信后端。
 *
 * 把串行实现单独放在一个文件里，使它与原生 MPI 后端的区别一目了然；集合操作
 * 保持相同接口，但退化为参数校验加一次直接拷贝。
 */
/*
 * Single-process communication backend.
 *
 * Keeping the serial implementation in its own file makes the distinction from
 * the native MPI backend explicit; collective operations preserve the same
 * interface but reduce to validation and a direct copy.
 */

#include "mpi_serial_backend.h"

#include "plasic_mpi.h"

#include <stddef.h>
#include <string.h>

/*
 * 串行启动：不调用任何 MPI 函数，直接报告“1 个进程、rank 0”。
 * 期望进程数不是 1 时返回 PLASIC_MP_START_SIZE_MISMATCH——这正是串行构建必须
 * 使用 NPRO=1 的原因。
 */
/*
 * Serial startup: calls no MPI function at all and simply reports "one
 * process, rank 0". An expected count other than 1 yields
 * PLASIC_MP_START_SIZE_MISMATCH, which is exactly why a serial build must use
 * NPRO=1.
 */
int32_t mpi_serial_start(
    int32_t expected_processes, int32_t *processes, int32_t *rank,
    int32_t *backend_error)
{
    if (processes == NULL || rank == NULL || backend_error == NULL)
    {
        return PLASIC_MP_START_ERROR;
    }
    *processes = 1;
    *rank = 0;
    *backend_error = 0;
    return expected_processes == 1 ? PLASIC_MP_START_OK
                                   : PLASIC_MP_START_SIZE_MISMATCH;
}

/* 单进程无需任何收尾工作。 */
/* Nothing needs to be torn down with one process. */
int32_t mpi_serial_stop(void)
{
    return 0;
}

/*
 * 广播 count 个 int32_t 元素：只有一个进程，数据本来就在“所有进程”手中，
 * 返回成功即可，buffer 保持不变。
 */
/*
 * Broadcast count int32_t elements: with one process the data is already in
 * every process's hands, so success is returned and buffer stays unchanged.
 */
int32_t mpi_serial_broadcast_integer(int32_t *buffer, int32_t count)
{
    (void)buffer;
    (void)count;
    return 0;
}

/* 广播 float 数组的串行版本，语义同上。 */
/* Serial version for a float array with the same semantics as above. */
int32_t mpi_serial_broadcast_real(float *buffer, int32_t count)
{
    (void)buffer;
    (void)count;
    return 0;
}

/*
 * 全局归约：单进程的“全局和”就是本地值，无需任何操作。
 */
/*
 * Global reduction: with one process the "global sum" is simply the local
 * value, so no operation is needed.
 */
int32_t mpi_serial_allreduce_real(float *buffer, int32_t count)
{
    (void)buffer;
    (void)count;
    return 0;
}

/*
 * 集合收集的串行退化：单进程时全局行数与局部行数必须相等（否则分块本身
 * 无意义），校验通过后一次 memcpy 完成“收集”。
 */
/*
 * Serial degradation of the collective gathers: with one process the global
 * and local row counts must be equal (otherwise the decomposition itself
 * would be meaningless); after that check a single memcpy performs the
 * "gather".
 */
static int32_t copy_collective(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    if (global_rows != local_rows)
    {
        return -1;
    }
    memcpy(global, local,
           (size_t)local_rows * (size_t)levels * sizeof(*global));
    return 0;
}

/*
 * 三个收集函数在单进程下语义完全相同，统一退化为 copy_collective。
 */
/*
 * All three gather functions have identical semantics with one process and
 * degrade to copy_collective.
 */
int32_t mpi_serial_allgather_grid(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    return copy_collective(
        global, local, global_rows, local_rows, levels);
}

int32_t mpi_serial_gather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    return copy_collective(
        global, local, global_rows, local_rows, levels);
}

int32_t mpi_serial_allgather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    return copy_collective(
        global, local, global_rows, local_rows, levels);
}

/*
 * 纬度 halo 交换：单个进程拥有全部纬度行，没有邻居需要交换；调用者预填的
 * 镜像边界值原样保留，因此串行与 MPI 构建的边界处理语义一致。
 */
/*
 * Latitude halo exchange: a single process owns every latitude row, so there
 * is no neighbour to exchange with; the caller's pre-filled mirrored edges are
 * kept unchanged, so the boundary handling has the same semantics in serial
 * and MPI builds.
 */
int32_t mpi_serial_exchange_latitude_halo(
    float *recv_north, float *recv_south, const float *send_north,
    const float *send_south, int32_t elements)
{
    (void)recv_north;
    (void)recv_south;
    (void)send_north;
    (void)send_south;
    (void)elements;
    return 0;
}
