#ifndef PLASIC_MPI_SERIAL_BACKEND_H
#define PLASIC_MPI_SERIAL_BACKEND_H

/*
 * 单进程后端的声明：不使用任何 MPI，把多进程接口按语义等价的方式退化为单进程
 * 操作，供 PLASIC_USE_MPI=0 的串行构建使用。
 *
 * 状态码约定与 mp_* 一致：成功返回 0，失败返回非零值。
 */
/*
 * Declarations of the single-process backend: it uses no MPI and degrades the
 * multi-process interface into semantically equivalent single-process
 * operations for serial builds with PLASIC_USE_MPI=0.
 *
 * Status convention matches mp_*: zero on success, non-zero on failure.
 */

#include <stdint.h>

/*
 * 串行启动：报告 1 个进程、rank 0；期望进程数不是 1 时返回不匹配状态。
 */
/*
 * Serial startup: reports one process and rank 0; returns a mismatch status
 * when the expected count is not 1.
 */
int32_t mpi_serial_start(
    int32_t expected_processes, int32_t *processes, int32_t *rank,
    int32_t *backend_error);
/* 单进程没有需要关闭的通信环境。 */
/* With one process there is no communication environment to shut down. */
int32_t mpi_serial_stop(void);

/* 广播与归约：单进程下数据本来就在本地，直接成功。 */
/* Broadcast and reduction: with one process the data is already local, so
 * they simply succeed. */
int32_t mpi_serial_broadcast_integer(int32_t *buffer, int32_t count);
int32_t mpi_serial_broadcast_real(float *buffer, int32_t count);
int32_t mpi_serial_allreduce_real(float *buffer, int32_t count);

/* 三个收集函数都退化为一次 memcpy。 */
/* All three gather functions degrade to a single memcpy. */
int32_t mpi_serial_allgather_grid(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);
int32_t mpi_serial_gather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);
int32_t mpi_serial_allgather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);

/* 纬度 halo 交换：单进程拥有全部纬度行，无需通信。 */
/* Latitude halo exchange: a single process owns every latitude row, so
 * nothing is communicated. */
int32_t mpi_serial_exchange_latitude_halo(
    float *recv_north, float *recv_south, const float *send_north,
    const float *send_south, int32_t elements);

#endif
