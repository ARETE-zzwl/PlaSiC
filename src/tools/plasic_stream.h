#ifndef PLASIC_STREAM_H
#define PLASIC_STREAM_H

#include <stdint.h>
#include <stdio.h>

enum plasic_stream_status
{
    PLASIC_STREAM_OK = 0,
    PLASIC_STREAM_BAD_ARGUMENT = 1,
    PLASIC_STREAM_IO_ERROR = 2,
    PLASIC_STREAM_UNSUPPORTED_PLATFORM = 3
};

enum plasic_stream_variable
{
    PLASIC_STREAM_AIR_TEMPERATURE = 1,
    PLASIC_STREAM_SPECIFIC_HUMIDITY = 2,
    PLASIC_STREAM_CLOUD_COVER = 3,
    PLASIC_STREAM_SURFACE_TEMPERATURE = 4,
    PLASIC_STREAM_PRECIPITATION = 5,
    PLASIC_STREAM_SEA_ICE_FRACTION = 6,
    PLASIC_STREAM_SEA_ICE_THICKNESS = 7,
    PLASIC_STREAM_SOIL_WATER = 8,
    PLASIC_STREAM_SEA_SURFACE_TEMPERATURE = 10,
    PLASIC_STREAM_ZONAL_WIND = 15,
    PLASIC_STREAM_MERIDIONAL_WIND = 16,
    PLASIC_STREAM_VERTICAL_VELOCITY = 17,
    PLASIC_STREAM_VORTICITY = 18,
    PLASIC_STREAM_DIVERGENCE = 19,
    PLASIC_STREAM_WIND_SPEED = 20,
    PLASIC_STREAM_SURFACE_PRESSURE = 21,
    PLASIC_STREAM_SURFACE_AIR_TEMPERATURE = 22,
    PLASIC_STREAM_SURFACE_HUMIDITY = 23,
    PLASIC_STREAM_EVAPORATION = 24,
    PLASIC_STREAM_SENSIBLE_HEAT_FLUX = 25,
    PLASIC_STREAM_LATENT_HEAT_FLUX = 26,
    PLASIC_STREAM_ALBEDO = 27,
    PLASIC_STREAM_SNOW_DEPTH = 28,
    PLASIC_STREAM_SOIL_TEMPERATURE = 30,
    PLASIC_STREAM_CONVECTIVE_PRECIPITATION = 31,
    PLASIC_STREAM_LARGE_SCALE_PRECIPITATION = 32,
    PLASIC_STREAM_SNOWFALL = 33,
    PLASIC_STREAM_ZONAL_WIND_STRESS = 34,
    PLASIC_STREAM_MERIDIONAL_WIND_STRESS = 35,
    PLASIC_STREAM_LAND_MASK = 36,
    PLASIC_STREAM_MIXED_LAYER_DEPTH = 37,
    PLASIC_STREAM_LIQUID_WATER = 38,
    PLASIC_STREAM_TOTAL_CLOUD_COVER = 39,
    PLASIC_STREAM_RELATIVE_HUMIDITY = 40,
    PLASIC_STREAM_GEOPOTENTIAL_HEIGHT = 41,
    PLASIC_STREAM_DEEP_OCEAN_TEMPERATURE = 42,
    PLASIC_STREAM_DEEP_OCEAN_SALINITY = 43
};

typedef struct plasic_stream_writer
{
    FILE *stream;
    int32_t nlon;
    int32_t nlat;
    int32_t nlev;
    int32_t ocean_levels;
} plasic_stream_writer;

int plasic_stream_open(
    plasic_stream_writer *writer, const char *path, int32_t nlon,
    int32_t nlat, int32_t nlev, int32_t ocean_levels);
int plasic_stream_write_frame(
    plasic_stream_writer *writer, int32_t variable_id, const char *name,
    const char *unit, int64_t step, const int32_t date_time[7],
    int32_t level_index, int32_t layers, const float *values);
int plasic_stream_close(plasic_stream_writer *writer);

#endif
