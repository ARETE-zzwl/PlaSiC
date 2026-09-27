#ifndef PLASIC_RESTARTMOD_IO_H
#define PLASIC_RESTARTMOD_IO_H

#include <stdint.h>

enum plasic_restart_status
{
    PLASIC_RESTART_OK = 0,
    PLASIC_RESTART_NOT_FOUND = 1,
    PLASIC_RESTART_IO_ERROR = 2,
    PLASIC_RESTART_BAD_RECORD = 3,
    PLASIC_RESTART_TOO_MANY_RECORDS = 4,
    PLASIC_RESTART_INVALID_ARGUMENT = 5,
    PLASIC_RESTART_COMMUNICATION_ERROR = 6,
    PLASIC_RESTART_ALLOCATION_ERROR = 7
};

/* Open and index a PlaSiC sequential-unformatted restart file. */
int32_t restart_open_read(const char *path);
int32_t restart_open_write(const char *path);
int32_t restart_close(void);

int32_t restart_read_integer(const char *name, int32_t *value);
int32_t restart_read_seed(const char *name, int32_t *values, int32_t count);
int32_t restart_read_array(
    const char *name, float *values, int32_t stored_rows, int32_t columns);
int32_t restart_read_array_slice(
    const char *name, float *values, int32_t global_rows, int32_t local_offset,
    int32_t local_rows, int32_t columns);

int32_t restart_write_seed(
    const char *name, const int32_t *values, int32_t count);

/* Coordinated restart writes: only root performs I/O and all ranks agree on
 * the returned status. Distributed inputs are gathered before writing. */
int32_t restart_write_root_integer(const char *name, int32_t value);
int32_t restart_write_root_array(
    const char *name, const float *values, int32_t stored_rows,
    int32_t columns);
int32_t restart_write_distributed_grid(
    const char *name, const float *local_values, int32_t global_rows,
    int32_t local_rows, int32_t columns);
int32_t restart_write_distributed_spectral(
    const char *name, const float *local_values, int32_t global_rows,
    int32_t local_rows, int32_t columns);

#endif
