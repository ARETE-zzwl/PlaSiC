#include "restartmod_io.h"
#include "mpi/plasic_mpi.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * PlaSiC's established restart format stores every logical record as
 *
 *     uint32 byte_count, byte_count payload bytes, uint32 byte_count
 *
 * in native byte order.  Each restart item consists of a 16-byte, blank-
 * padded name record followed by one data record.  The standalone C runtime
 * owns both binary I/O and the read-side record index.
 */

enum
{
    RESTART_NAME_BYTES = 16,
    RESTART_MAX_RECORDS = 200
};

typedef struct restart_record
{
    char name[RESTART_NAME_BYTES];
    long data_offset;
    uint32_t data_bytes;
} restart_record;

static FILE *read_file;
static FILE *write_file;
static restart_record records[RESTART_MAX_RECORDS];
static int32_t record_count;

static int32_t broadcast_restart_status(int32_t status)
{
    return mp_broadcast_integer(&status, 1) == 0
               ? status
               : PLASIC_RESTART_COMMUNICATION_ERROR;
}

_Static_assert(sizeof(int32_t) == 4, "restart INTEGER records require 32 bits");
_Static_assert(sizeof(float) == 4, "restart REAL records require 32 bits");
_Static_assert(sizeof(uint32_t) == 4, "restart markers require 32 bits");

static char *copy_path(const char *path)
{
    size_t path_length;
    char *copy;

    if (path == NULL || path[0] == '\0')
    {
        return NULL;
    }
    path_length = strlen(path);
    copy = (char *)malloc(path_length + 1);
    if (copy == NULL)
    {
        return NULL;
    }
    memcpy(copy, path, path_length + 1);
    return copy;
}

static void make_name(char result[RESTART_NAME_BYTES], const char *name)
{
    size_t copied = strlen(name);

    memset(result, ' ', RESTART_NAME_BYTES);
    if (copied > RESTART_NAME_BYTES)
    {
        copied = RESTART_NAME_BYTES;
    }
    if (copied != 0)
    {
        memcpy(result, name, copied);
    }
}

static int read_marker(FILE *stream, uint32_t *marker, int *at_eof)
{
    size_t count;

    *at_eof = 0;
    count = fread(marker, 1, sizeof(*marker), stream);
    if (count == sizeof(*marker))
    {
        return PLASIC_RESTART_OK;
    }
    if (count == 0 && feof(stream))
    {
        *at_eof = 1;
        return PLASIC_RESTART_OK;
    }
    return PLASIC_RESTART_IO_ERROR;
}

static int write_marker(FILE *stream, uint32_t marker)
{
    if (stream == NULL || fwrite(&marker, sizeof(marker), 1, stream) != 1)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    return PLASIC_RESTART_OK;
}

static int close_stream(FILE **stream)
{
    int result = PLASIC_RESTART_OK;

    if (*stream != NULL)
    {
        if (fclose(*stream) != 0)
        {
            result = PLASIC_RESTART_IO_ERROR;
        }
        *stream = NULL;
    }
    return result;
}

static int scan_records(void)
{
    uint32_t leading;
    uint32_t trailing;
    int at_eof;
    int status;

    record_count = 0;
    for (;;)
    {
        char name[RESTART_NAME_BYTES];
        long data_offset;

        status = read_marker(read_file, &leading, &at_eof);
        if (status != PLASIC_RESTART_OK)
        {
            return status;
        }
        if (at_eof)
        {
            return PLASIC_RESTART_OK;
        }
        if (leading != RESTART_NAME_BYTES)
        {
            return PLASIC_RESTART_BAD_RECORD;
        }
        if (fread(name, 1, RESTART_NAME_BYTES, read_file) !=
            RESTART_NAME_BYTES)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        status = read_marker(read_file, &trailing, &at_eof);
        if (status != PLASIC_RESTART_OK || at_eof)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        if (trailing != leading)
        {
            return PLASIC_RESTART_BAD_RECORD;
        }

        status = read_marker(read_file, &leading, &at_eof);
        if (status != PLASIC_RESTART_OK || at_eof)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        data_offset = ftell(read_file);
        if (data_offset < 0)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        if (leading > (uint32_t)LONG_MAX ||
            fseek(read_file, (long)leading, SEEK_CUR) != 0)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        status = read_marker(read_file, &trailing, &at_eof);
        if (status != PLASIC_RESTART_OK || at_eof)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        if (trailing != leading)
        {
            return PLASIC_RESTART_BAD_RECORD;
        }

        memcpy(records[record_count].name, name, RESTART_NAME_BYTES);
        records[record_count].data_offset = data_offset;
        records[record_count].data_bytes = leading;
        ++record_count;

        /* The original code stops as soon as record number 200 is seen. */
        if (record_count >= RESTART_MAX_RECORDS)
        {
            return PLASIC_RESTART_TOO_MANY_RECORDS;
        }
    }
}

static restart_record *find_record(const char *name)
{
    char padded[RESTART_NAME_BYTES];
    int32_t index;

    if (name == NULL)
    {
        return NULL;
    }
    make_name(padded, name);
    for (index = 0; index < record_count; ++index)
    {
        if (memcmp(records[index].name, padded, RESTART_NAME_BYTES) == 0)
        {
            return &records[index];
        }
    }
    return NULL;
}

static int write_name_record(const char *name)
{
    char padded[RESTART_NAME_BYTES];

    if (name == NULL)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    make_name(padded, name);
    if (write_marker(write_file, RESTART_NAME_BYTES) != PLASIC_RESTART_OK)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    if (fwrite(padded, 1, RESTART_NAME_BYTES, write_file) !=
        RESTART_NAME_BYTES)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    return write_marker(write_file, RESTART_NAME_BYTES);
}

static int read_contiguous(const char *name,
                           void *destination, uint32_t required_bytes)
{
    restart_record *record;

    if (read_file == NULL || (destination == NULL && required_bytes != 0))
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    record = find_record(name);
    if (record == NULL)
    {
        return PLASIC_RESTART_NOT_FOUND;
    }
    if (record->data_bytes < required_bytes)
    {
        return PLASIC_RESTART_BAD_RECORD;
    }
    if (fseek(read_file, record->data_offset, SEEK_SET) != 0)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    if (required_bytes != 0 &&
        fread(destination, 1, required_bytes, read_file) != required_bytes)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    return PLASIC_RESTART_OK;
}

static int checked_payload_bytes(int32_t count, size_t element_bytes,
                                 uint32_t *payload_bytes)
{
    size_t total;

    if (count < 0)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    if ((size_t)count > UINT32_MAX / element_bytes)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    total = (size_t)count * element_bytes;
    *payload_bytes = (uint32_t)total;
    return PLASIC_RESTART_OK;
}

int32_t restart_open_read(const char *path)
{
    char *path_copy;
    int status;

    path_copy = copy_path(path);
    if (path_copy == NULL)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    (void)close_stream(&read_file);
    read_file = fopen(path_copy, "rb");
    free(path_copy);
    if (read_file == NULL)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    status = scan_records();
    if (status != PLASIC_RESTART_OK &&
        status != PLASIC_RESTART_TOO_MANY_RECORDS)
    {
        (void)close_stream(&read_file);
    }
    return status;
}

int32_t restart_open_write(const char *path)
{
    char *path_copy;

    path_copy = copy_path(path);
    if (path_copy == NULL)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    (void)close_stream(&write_file);
    write_file = fopen(path_copy, "wb");
    free(path_copy);
    if (write_file == NULL)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    return PLASIC_RESTART_OK;
}

int32_t restart_close(void)
{
    int read_status = close_stream(&read_file);
    int write_status = close_stream(&write_file);

    record_count = 0;
    if (read_status != PLASIC_RESTART_OK)
    {
        return read_status;
    }
    return write_status;
}

int32_t restart_read_integer(const char *name, int32_t *value)
{
    return read_contiguous(name, value, sizeof(*value));
}

int32_t restart_read_seed(const char *name, int32_t *values, int32_t count)
{
    uint32_t payload_bytes;
    int status = checked_payload_bytes(count, sizeof(*values), &payload_bytes);

    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    return read_contiguous(name, values, payload_bytes);
}

int32_t restart_read_array(const char *name, float *values,
                           int32_t stored_rows, int32_t columns)
{
    restart_record *record;
    uint32_t payload_bytes;
    int32_t column;
    int status;

    if (stored_rows < 0 || columns < 0 ||
        (values == NULL && stored_rows != 0 && columns != 0))
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    if (columns != 0 && stored_rows > INT32_MAX / columns)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    status = checked_payload_bytes(stored_rows * columns, sizeof(*values),
                                   &payload_bytes);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    record = find_record(name);
    if (read_file == NULL || record == NULL)
    {
        return record == NULL ? PLASIC_RESTART_NOT_FOUND
                              : PLASIC_RESTART_IO_ERROR;
    }
    if (record->data_bytes < payload_bytes)
    {
        return PLASIC_RESTART_BAD_RECORD;
    }
    if (fseek(read_file, record->data_offset, SEEK_SET) != 0)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    if (payload_bytes == 0)
    {
        return PLASIC_RESTART_OK;
    }
    for (column = 0; column < columns; ++column)
    {
        float *column_start = values +
            (size_t)column * (size_t)stored_rows;
        if (stored_rows != 0 &&
            fread(column_start, sizeof(*values), (size_t)stored_rows,
                  read_file) != (size_t)stored_rows)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
    }
    return PLASIC_RESTART_OK;
}

int32_t restart_read_array_slice(
    const char *name, float *values,
    int32_t global_rows, int32_t local_offset, int32_t local_rows,
    int32_t columns)
{
    restart_record *record;
    uint32_t expected_bytes;
    int32_t column;
    int32_t rows_to_read;
    int status;

    if (global_rows < 0 || local_offset < 0 || local_rows < 0 ||
        columns < 0 ||
        local_offset > global_rows ||
        (values == NULL && local_rows != 0 && columns != 0))
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    if (columns != 0 && global_rows > INT32_MAX / columns)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    status = checked_payload_bytes(
        global_rows * columns, sizeof(*values), &expected_bytes);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    record = find_record(name);
    if (read_file == NULL || record == NULL)
    {
        return record == NULL ? PLASIC_RESTART_NOT_FOUND
                              : PLASIC_RESTART_IO_ERROR;
    }
    if (record->data_bytes < expected_bytes)
    {
        return PLASIC_RESTART_BAD_RECORD;
    }
    rows_to_read = local_rows;
    if (rows_to_read > global_rows - local_offset)
    {
        rows_to_read = global_rows - local_offset;
    }

    for (column = 0; column < columns; ++column)
    {
        float *column_start =
            values + (size_t)column * (size_t)local_rows;
        const long byte_offset =
            (long)(((size_t)column * (size_t)global_rows +
                    (size_t)local_offset) * sizeof(*values));

        if (fseek(
                read_file, record->data_offset + byte_offset,
                SEEK_SET) != 0)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        if (rows_to_read != 0 &&
            fread(column_start, sizeof(*values), (size_t)rows_to_read,
                  read_file) != (size_t)rows_to_read)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
        if (rows_to_read < local_rows)
        {
            memset(
                column_start + rows_to_read, 0,
                (size_t)(local_rows - rows_to_read) * sizeof(*values));
        }
    }
    return PLASIC_RESTART_OK;
}

static int32_t restart_write_integer(const char *name, int32_t value)
{
    const uint32_t payload_bytes = sizeof(value);
    int status = write_name_record(name);

    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    status = write_marker(write_file, payload_bytes);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    if (fwrite(&value, sizeof(value), 1, write_file) != 1)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    return write_marker(write_file, payload_bytes);
}
int32_t restart_write_seed(const char *name,
                           const int32_t *values, int32_t count)
{
    uint32_t payload_bytes;
    int status;

    if (values == NULL && count != 0)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    status = checked_payload_bytes(count, sizeof(*values), &payload_bytes);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    status = write_name_record(name);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    status = write_marker(write_file, payload_bytes);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    if (count != 0 &&
        fwrite(values, sizeof(*values), (size_t)count, write_file) !=
            (size_t)count)
    {
        return PLASIC_RESTART_IO_ERROR;
    }
    return write_marker(write_file, payload_bytes);
}

static int32_t restart_write_array(const char *name, const float *values,
                                   int32_t stored_rows, int32_t columns)
{
    uint32_t payload_bytes;
    int32_t column;
    int status;

    if (stored_rows < 0 || columns < 0 ||
        (values == NULL && stored_rows != 0 && columns != 0))
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    if (columns != 0 && stored_rows > INT32_MAX / columns)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    status = checked_payload_bytes(stored_rows * columns, sizeof(*values),
                                   &payload_bytes);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    status = write_name_record(name);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    status = write_marker(write_file, payload_bytes);
    if (status != PLASIC_RESTART_OK)
    {
        return status;
    }
    if (payload_bytes == 0)
    {
        return write_marker(write_file, payload_bytes);
    }
    for (column = 0; column < columns; ++column)
    {
        const float *column_start =
            values + (size_t)column * (size_t)stored_rows;
        if (stored_rows != 0 &&
            fwrite(column_start, sizeof(*values), (size_t)stored_rows,
                   write_file) != (size_t)stored_rows)
        {
            return PLASIC_RESTART_IO_ERROR;
        }
    }
    return write_marker(write_file, payload_bytes);
}

int32_t restart_write_root_integer(const char *name, int32_t value)
{
    int32_t status = PLASIC_RESTART_OK;

    if (name == NULL)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    if (mp_is_root())
    {
        status = restart_write_integer(name, value);
    }
    return broadcast_restart_status(status);
}

int32_t restart_write_root_array(
    const char *name, const float *values, int32_t stored_rows,
    int32_t columns)
{
    int32_t status = PLASIC_RESTART_OK;

    if (name == NULL)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    if (mp_is_root())
    {
        status = restart_write_array(name, values, stored_rows, columns);
    }
    return broadcast_restart_status(status);
}

static float *allocate_gather_buffer(int32_t rows, int32_t columns)
{
    size_t values;

    if (rows < 0 || columns < 0 ||
        (rows != 0 && (size_t)columns > SIZE_MAX / (size_t)rows))
    {
        return NULL;
    }
    values = (size_t)rows * (size_t)columns;
    if (values > SIZE_MAX / sizeof(float))
    {
        return NULL;
    }
    return (float *)malloc((values == 0U ? 1U : values) * sizeof(float));
}

int32_t restart_write_distributed_grid(
    const char *name, const float *local_values, int32_t global_rows,
    int32_t local_rows, int32_t columns)
{
    float *global_values;
    int32_t status;

    if (name == NULL || local_values == NULL || global_rows < 0 ||
        local_rows < 0 || columns < 0)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    global_values = allocate_gather_buffer(global_rows, columns);
    if (global_values == NULL)
    {
        return PLASIC_RESTART_ALLOCATION_ERROR;
    }

    status = mp_allgather_grid(
                 global_values, local_values, global_rows, local_rows,
                 columns) == 0
                 ? PLASIC_RESTART_OK
                 : PLASIC_RESTART_COMMUNICATION_ERROR;
    if (status == PLASIC_RESTART_OK)
    {
        status = restart_write_root_array(
            name, global_values, global_rows, columns);
    }
    free(global_values);
    return status;
}

int32_t restart_write_distributed_spectral(
    const char *name, const float *local_values, int32_t global_rows,
    int32_t local_rows, int32_t columns)
{
    float *global_values;
    int32_t status;

    if (name == NULL || local_values == NULL || global_rows < 0 ||
        local_rows < 0 || columns < 0)
    {
        return PLASIC_RESTART_INVALID_ARGUMENT;
    }
    global_values = allocate_gather_buffer(global_rows, columns);
    if (global_values == NULL)
    {
        return PLASIC_RESTART_ALLOCATION_ERROR;
    }

    status = mp_gather_spectral(
                 global_values, local_values, global_rows,
                 local_rows, columns) == 0
                 ? PLASIC_RESTART_OK
                 : PLASIC_RESTART_COMMUNICATION_ERROR;
    if (status == PLASIC_RESTART_OK)
    {
        status = restart_write_root_array(
            name, global_values, global_rows, columns);
    }
    free(global_values);
    return status;
}
