#ifndef PLASIC_SURFACE_DATA_IO_H
#define PLASIC_SURFACE_DATA_IO_H

#include <stdint.h>

enum surface_data_status
{
    SURFACE_DATA_OK = 0,
    SURFACE_DATA_NOT_FOUND = 1,
    SURFACE_DATA_IO_ERROR = 2,
    SURFACE_DATA_BAD_HEADER = 3,
    SURFACE_DATA_BAD_ARGUMENT = 4,
    SURFACE_DATA_WRONG_MONTH_COUNT = 5,
    SURFACE_DATA_FIELD_MISMATCH = 6,
    SURFACE_DATA_GRID_MISMATCH = 7
};

/*
 * Read one PlaSiC Nxxx_surf_<field_name>.nc field and retain the current
 * process's contiguous, complete-latitude slice. Runtime input is NetCDF
 * only: there is intentionally no SRA fallback. A one-record field is
 * repeated when 14 records are requested; a 12-month field receives cyclic
 * December/January halos, while an already expanded 14-record field is
 * retained verbatim.
 */
int surface_data_read_local(
    const char *data_directory, int32_t truncation, int32_t latitude_count,
    int32_t longitude_count, const char *field_name, int32_t requested_months,
    int32_t local_offset, int32_t local_count, float *values);

#endif
