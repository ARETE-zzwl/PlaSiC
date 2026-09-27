#include "surface_data_io.h"

#include "surfmod_kernels.h"

#include <errno.h>
#include <limits.h>
#include <netcdf.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

enum
{
    SURFACE_PATH_BYTES = 1024,
    SURFACE_FIELD_NAME_BYTES = 64,
    SURFACE_VARIABLE_RANK = 3
};

static int field_name_is_valid(const char *field_name)
{
    if (field_name == NULL || field_name[0] == '\0' ||
        strlen(field_name) >= SURFACE_FIELD_NAME_BYTES)
    {
        return 0;
    }
    return strchr(field_name, '/') == NULL;
}

static int build_path(
    char path[SURFACE_PATH_BYTES], const char *data_directory,
    int32_t truncation, int32_t latitude_count, const char *field_name)
{
    const int count = snprintf(
        path, SURFACE_PATH_BYTES, "%s/T%d/N%03d_surf_%s.nc",
        data_directory, truncation, latitude_count, field_name);

    return count >= 0 && count < SURFACE_PATH_BYTES;
}

static int open_status(int netcdf_status)
{
    /* NetCDF may return the platform ENOENT directly for nc_open(). */
    if (netcdf_status == ENOENT)
    {
        return SURFACE_DATA_NOT_FOUND;
    }
#ifdef NC_ENOTFOUND
    if (netcdf_status == NC_ENOTFOUND)
    {
        return SURFACE_DATA_NOT_FOUND;
    }
#endif
    return SURFACE_DATA_IO_ERROR;
}

static int read_required_global_integer(
    int file_id, const char *name, int32_t *value)
{
    int result;

    if (nc_get_att_int(file_id, NC_GLOBAL, name, &result) != NC_NOERR)
    {
        return 0;
    }
    *value = (int32_t)result;
    return 1;
}

static int read_required_attribute(
    int file_id, int variable_id, const char *name,
    char value[SURFACE_FIELD_NAME_BYTES])
{
    size_t length;

    if (nc_inq_attlen(file_id, variable_id, name, &length) != NC_NOERR ||
        length == 0 || length >= SURFACE_FIELD_NAME_BYTES ||
        nc_get_att_text(file_id, variable_id, name, value) != NC_NOERR)
    {
        return 0;
    }
    value[length] = '\0';
    return 1;
}

int surface_data_read_local(
    const char *data_directory, int32_t truncation, int32_t latitude_count,
    int32_t longitude_count, const char *field_name, int32_t requested_months,
    int32_t local_offset, int32_t local_count, float *values)
{
    char path[SURFACE_PATH_BYTES];
    int file_id = -1;
    int variable_id;
    int dimension_ids[SURFACE_VARIABLE_RANK];
    int expected_dimension_ids[SURFACE_VARIABLE_RANK];
    int variable_dimensions;
    char stored_name[SURFACE_FIELD_NAME_BYTES];
    char attribute_name[SURFACE_FIELD_NAME_BYTES];
    int dimension;
    int32_t stored_truncation;
    int32_t stored_latitudes;
    int32_t stored_longitudes;
    size_t dimension_lengths[SURFACE_VARIABLE_RANK];
    size_t start[SURFACE_VARIABLE_RANK];
    size_t count[SURFACE_VARIABLE_RANK];
    int32_t source_months;
    int32_t months_to_read;
    int status = SURFACE_DATA_OK;
    int netcdf_status;

    if (data_directory == NULL || values == NULL ||
        !field_name_is_valid(field_name) || truncation <= 0 ||
        latitude_count <= 0 || longitude_count <= 0 ||
        (requested_months != 1 && requested_months != 14) ||
        local_offset < 0 || local_count < 0 ||
        latitude_count > INT32_MAX / longitude_count ||
        local_offset % longitude_count != 0 ||
        local_count % longitude_count != 0)
    {
        return SURFACE_DATA_BAD_ARGUMENT;
    }
    if (local_offset > latitude_count * longitude_count ||
        local_count > latitude_count * longitude_count - local_offset)
    {
        return SURFACE_DATA_BAD_ARGUMENT;
    }
    if (!build_path(path, data_directory, truncation, latitude_count, field_name))
    {
        return SURFACE_DATA_BAD_ARGUMENT;
    }

    netcdf_status = nc_open(path, NC_NOWRITE, &file_id);
    if (netcdf_status != NC_NOERR)
    {
        return open_status(netcdf_status);
    }

    if (!read_required_attribute(file_id, NC_GLOBAL, "field_name", stored_name) ||
        nc_inq_varid(file_id, "surface", &variable_id) != NC_NOERR ||
        !read_required_attribute(
            file_id, variable_id, "field_name", attribute_name))
    {
        status = SURFACE_DATA_BAD_HEADER;
        goto done;
    }
    if (strcmp(stored_name, field_name) != 0 ||
        strcmp(attribute_name, field_name) != 0)
    {
        status = SURFACE_DATA_FIELD_MISMATCH;
        goto done;
    }
    if (!read_required_global_integer(
            file_id, "spectral_truncation", &stored_truncation) ||
        !read_required_global_integer(
            file_id, "latitude_count", &stored_latitudes) ||
        !read_required_global_integer(
            file_id, "longitude_count", &stored_longitudes))
    {
        status = SURFACE_DATA_BAD_HEADER;
        goto done;
    }
    if (stored_truncation != truncation ||
        stored_latitudes != latitude_count ||
        stored_longitudes != longitude_count)
    {
        status = SURFACE_DATA_GRID_MISMATCH;
        goto done;
    }

    if (nc_inq_dimid(file_id, "time", &expected_dimension_ids[0]) != NC_NOERR ||
        nc_inq_dimid(file_id, "lat", &expected_dimension_ids[1]) != NC_NOERR ||
        nc_inq_dimid(file_id, "lon", &expected_dimension_ids[2]) != NC_NOERR ||
        nc_inq_varndims(file_id, variable_id, &variable_dimensions) != NC_NOERR ||
        variable_dimensions != SURFACE_VARIABLE_RANK ||
        nc_inq_vardimid(file_id, variable_id, dimension_ids) != NC_NOERR)
    {
        status = SURFACE_DATA_BAD_HEADER;
        goto done;
    }
    for (dimension = 0; dimension < SURFACE_VARIABLE_RANK; ++dimension)
    {
        if (dimension_ids[dimension] != expected_dimension_ids[dimension] ||
            nc_inq_dimlen(
                file_id, dimension_ids[dimension],
                &dimension_lengths[dimension]) != NC_NOERR)
        {
            status = SURFACE_DATA_BAD_HEADER;
            goto done;
        }
    }
    if (dimension_lengths[1] != (size_t)latitude_count ||
        dimension_lengths[2] != (size_t)longitude_count)
    {
        status = SURFACE_DATA_GRID_MISMATCH;
        goto done;
    }
    if (dimension_lengths[0] > (size_t)INT32_MAX)
    {
        status = SURFACE_DATA_WRONG_MONTH_COUNT;
        goto done;
    }
    source_months = (int32_t)dimension_lengths[0];
    if (source_months <= 0 ||
        (requested_months == 14 && source_months != 1 &&
         source_months != 12 && source_months != 14))
    {
        status = SURFACE_DATA_WRONG_MONTH_COUNT;
        goto done;
    }
    months_to_read = requested_months == 1 ? 1 : source_months;
    start[0] = 0;
    start[1] = (size_t)(local_offset / longitude_count);
    start[2] = 0;
    count[0] = (size_t)months_to_read;
    count[1] = (size_t)(local_count / longitude_count);
    count[2] = (size_t)longitude_count;
    if (nc_get_vara_float(file_id, variable_id, start, count, values) !=
        NC_NOERR)
    {
        status = SURFACE_DATA_IO_ERROR;
        goto done;
    }
    if (requested_months == 14 && source_months != 14)
    {
        surface_expand_months(local_count, source_months, values);
    }

done:
    if (nc_close(file_id) != NC_NOERR && status == SURFACE_DATA_OK)
    {
        status = SURFACE_DATA_IO_ERROR;
    }
    return status;
}
