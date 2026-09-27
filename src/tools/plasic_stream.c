#include "plasic_stream.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

enum
{
    PLASIC_STREAM_VERSION = 1,
    PLASIC_STREAM_FILE_HEADER_BYTES = 64,
    PLASIC_STREAM_FRAME_HEADER_BYTES = 128
};

static void store_u16(unsigned char *buffer, size_t offset, uint16_t value)
{
    buffer[offset] = (unsigned char)(value & UINT16_C(0xff));
    buffer[offset + 1] =
        (unsigned char)((value >> 8) & UINT16_C(0xff));
}

static void store_u32(unsigned char *buffer, size_t offset, uint32_t value)
{
    size_t byte;

    for (byte = 0; byte < 4; ++byte)
    {
        buffer[offset + byte] =
            (unsigned char)((value >> (byte * 8)) & UINT32_C(0xff));
    }
}

static void store_u64(unsigned char *buffer, size_t offset, uint64_t value)
{
    size_t byte;

    for (byte = 0; byte < 8; ++byte)
    {
        buffer[offset + byte] =
            (unsigned char)((value >> (byte * 8)) & UINT64_C(0xff));
    }
}

static void store_float(
    unsigned char *buffer, size_t offset, float value)
{
    uint32_t bits;

    memcpy(&bits, &value, sizeof(bits));
    store_u32(buffer, offset, bits);
}

static uint32_t fnv1a(const unsigned char *bytes, size_t count)
{
    uint32_t hash = UINT32_C(2166136261);
    size_t index;

    for (index = 0; index < count; ++index)
    {
        hash ^= bytes[index];
        hash *= UINT32_C(16777619);
    }
    return hash;
}

static int platform_is_supported(void)
{
    const uint16_t marker = UINT16_C(1);
    const unsigned char *bytes = (const unsigned char *)&marker;

    return sizeof(float) == 4 && bytes[0] == 1;
}

int plasic_stream_open(
    plasic_stream_writer *writer, const char *path, int32_t nlon,
    int32_t nlat, int32_t nlev, int32_t ocean_levels)
{
    unsigned char header[PLASIC_STREAM_FILE_HEADER_BYTES] = {0};

    if (writer == NULL || path == NULL || path[0] == '\0' ||
        nlon <= 0 || nlat <= 0 || nlev <= 0 || ocean_levels <= 0)
    {
        return PLASIC_STREAM_BAD_ARGUMENT;
    }
    if (!platform_is_supported())
    {
        return PLASIC_STREAM_UNSUPPORTED_PLATFORM;
    }
    memset(writer, 0, sizeof(*writer));
    writer->stream = fopen(path, "wb");
    if (writer->stream == NULL)
    {
        return PLASIC_STREAM_IO_ERROR;
    }
    memcpy(header, "PLASICPS", 8);
    store_u16(header, 8, PLASIC_STREAM_VERSION);
    store_u16(header, 10, PLASIC_STREAM_FILE_HEADER_BYTES);
    store_u32(header, 12, UINT32_C(0));
    store_u32(header, 16, UINT32_C(0x01020304));
    store_u32(header, 20, (uint32_t)nlon);
    store_u32(header, 24, (uint32_t)nlat);
    store_u32(header, 28, (uint32_t)nlev);
    store_u32(header, 32, (uint32_t)ocean_levels);
    store_u32(header, 36, (uint32_t)sizeof(float));
    if (fwrite(header, 1, sizeof(header), writer->stream) !=
        sizeof(header))
    {
        (void)fclose(writer->stream);
        writer->stream = NULL;
        return PLASIC_STREAM_IO_ERROR;
    }
    writer->nlon = nlon;
    writer->nlat = nlat;
    writer->nlev = nlev;
    writer->ocean_levels = ocean_levels;
    return PLASIC_STREAM_OK;
}

int plasic_stream_write_frame(
    plasic_stream_writer *writer, int32_t variable_id, const char *name,
    const char *unit, int64_t step, const int32_t date_time[7],
    int32_t level_index, int32_t layers, const float *values)
{
    unsigned char header[PLASIC_STREAM_FRAME_HEADER_BYTES] = {0};
    size_t value_count;
    size_t payload_bytes;
    float minimum;
    float maximum;
    size_t index;

    if (writer == NULL || writer->stream == NULL || variable_id <= 0 ||
        name == NULL || unit == NULL || date_time == NULL ||
        layers <= 0 || values == NULL ||
        (size_t)writer->nlon > SIZE_MAX / (size_t)writer->nlat)
    {
        return PLASIC_STREAM_BAD_ARGUMENT;
    }
    value_count =
        (size_t)writer->nlon * (size_t)writer->nlat;
    if (value_count > SIZE_MAX / (size_t)layers)
    {
        return PLASIC_STREAM_BAD_ARGUMENT;
    }
    value_count *= (size_t)layers;
    if (value_count > SIZE_MAX / sizeof(*values))
    {
        return PLASIC_STREAM_BAD_ARGUMENT;
    }
    payload_bytes = value_count * sizeof(*values);
    if (payload_bytes > UINT32_MAX)
    {
        return PLASIC_STREAM_BAD_ARGUMENT;
    }
    minimum = values[0];
    maximum = values[0];
    for (index = 0; index < value_count; ++index)
    {
        if (!isfinite(values[index]))
        {
            return PLASIC_STREAM_BAD_ARGUMENT;
        }
        minimum = fminf(minimum, values[index]);
        maximum = fmaxf(maximum, values[index]);
    }

    memcpy(header, "PFRAME1\0", 8);
    store_u16(header, 8, PLASIC_STREAM_VERSION);
    store_u16(header, 10, PLASIC_STREAM_FRAME_HEADER_BYTES);
    store_u32(header, 12, (uint32_t)variable_id);
    store_u32(header, 16, (uint32_t)payload_bytes);
    store_u64(header, 20, (uint64_t)step);
    store_u32(header, 28, (uint32_t)date_time[0]);
    store_u32(header, 32, (uint32_t)date_time[1]);
    store_u32(header, 36, (uint32_t)date_time[2]);
    store_u32(header, 40, (uint32_t)date_time[3]);
    store_u32(header, 44, (uint32_t)date_time[4]);
    store_u32(header, 48, (uint32_t)level_index);
    store_u32(header, 52, (uint32_t)writer->nlon);
    store_u32(header, 56, (uint32_t)writer->nlat);
    store_u32(header, 60, (uint32_t)layers);
    store_float(header, 64, minimum);
    store_float(header, 68, maximum);
    store_u32(
        header, 72,
        fnv1a((const unsigned char *)values, payload_bytes));
    (void)snprintf((char *)header + 76, 32, "%s", name);
    (void)snprintf((char *)header + 108, 16, "%s", unit);

    if (fwrite(header, 1, sizeof(header), writer->stream) !=
            sizeof(header) ||
        fwrite(values, 1, payload_bytes, writer->stream) !=
            payload_bytes ||
        fflush(writer->stream) != 0)
    {
        return PLASIC_STREAM_IO_ERROR;
    }
    return PLASIC_STREAM_OK;
}

int plasic_stream_close(plasic_stream_writer *writer)
{
    int status = PLASIC_STREAM_OK;

    if (writer == NULL)
    {
        return PLASIC_STREAM_BAD_ARGUMENT;
    }
    if (writer->stream != NULL && fclose(writer->stream) != 0)
    {
        status = PLASIC_STREAM_IO_ERROR;
    }
    memset(writer, 0, sizeof(*writer));
    return status;
}
