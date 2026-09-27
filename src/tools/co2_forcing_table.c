#include "co2_forcing_table.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum
{
    CO2_FORCING_LINE_BYTES = 512
};

static void clear_error(co2_forcing_table_error *error)
{
    if (error != NULL)
    {
        error->year = 0;
        error->month = 0;
        error->band = 0;
        error->latitude = 0.0f;
        error->first_year = 0;
        error->last_year = 0;
    }
}

static int parse_record(
    const char *line, float *value, int32_t *year,
    int32_t *month, float *latitude)
{
    int parsed_year;
    int parsed_month;

    if (sscanf(
            line, "%*[^,],%*[^,],%f,%d,%d,%f",
            value, &parsed_year, &parsed_month, latitude) != 4 ||
        parsed_year < 0 || parsed_year > INT32_MAX ||
        parsed_month < 1 || parsed_month > CO2_FORCING_MONTHS ||
        !isfinite(*value) || *value <= 0.0f ||
        !isfinite(*latitude))
    {
        return 0;
    }
    *year = (int32_t)parsed_year;
    *month = (int32_t)parsed_month;
    return 1;
}

static int latitude_band(float latitude)
{
    const float fractional = (latitude + 82.5f) / 15.0f;
    const int band = (int)lroundf(fractional);
    const float expected = -82.5f + 15.0f * (float)band;

    return band >= 0 && band < CO2_FORCING_LATITUDE_BANDS &&
                   fabsf(latitude - expected) < 1.0e-3f
               ? band
               : -1;
}

int co2_forcing_table_load(
    const char *path, co2_forcing_table *table,
    co2_forcing_table_error *error)
{
    unsigned char *seen = NULL;
    float *values = NULL;
    char line[CO2_FORCING_LINE_BYTES];
    FILE *stream = NULL;
    int32_t first_year = INT32_MAX;
    int32_t last_year = -1;
    size_t value_count;

    clear_error(error);
    table->values = NULL;
    table->first_year = 0;
    table->last_year = 0;

    stream = fopen(path, "r");
    if (stream == NULL)
    {
        return CO2_FORCING_TABLE_OPEN_FAILED;
    }
    if (fgets(line, sizeof(line), stream) == NULL ||
        strstr(line, "data") == NULL || strstr(line, "year") == NULL ||
        strstr(line, "month") == NULL || strstr(line, "lat") == NULL)
    {
        (void)fclose(stream);
        return CO2_FORCING_TABLE_BAD_HEADER;
    }
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        float value;
        float latitude;
        int32_t year;
        int32_t month;

        if (!parse_record(line, &value, &year, &month, &latitude) ||
            latitude_band(latitude) < 0)
        {
            (void)fclose(stream);
            return CO2_FORCING_TABLE_BAD_RECORD;
        }
        first_year = year < first_year ? year : first_year;
        last_year = year > last_year ? year : last_year;
    }
    if (ferror(stream) || first_year == INT32_MAX || last_year < first_year)
    {
        (void)fclose(stream);
        return CO2_FORCING_TABLE_SCAN_FAILED;
    }
    value_count = ((size_t)(last_year - first_year) + 1) *
                  CO2_FORCING_MONTHS * CO2_FORCING_LATITUDE_BANDS;
    if (value_count > (size_t)INT32_MAX)
    {
        (void)fclose(stream);
        return CO2_FORCING_TABLE_TOO_LARGE;
    }
    values = (float *)malloc(value_count * sizeof(*values));
    seen = (unsigned char *)calloc(value_count, sizeof(*seen));
    if (values == NULL || seen == NULL)
    {
        free(values);
        free(seen);
        (void)fclose(stream);
        return CO2_FORCING_TABLE_ALLOCATION_FAILED;
    }
    if (fseek(stream, 0, SEEK_SET) != 0 ||
        fgets(line, sizeof(line), stream) == NULL)
    {
        free(values);
        free(seen);
        (void)fclose(stream);
        return CO2_FORCING_TABLE_REWIND_FAILED;
    }
    while (fgets(line, sizeof(line), stream) != NULL)
    {
        float value;
        float latitude;
        int32_t year;
        int32_t month;
        int band;
        size_t index;

        if (!parse_record(line, &value, &year, &month, &latitude))
        {
            free(values);
            free(seen);
            (void)fclose(stream);
            return CO2_FORCING_TABLE_BAD_RECORD;
        }
        band = latitude_band(latitude);
        if (band < 0)
        {
            if (error != NULL)
            {
                error->latitude = latitude;
            }
            free(values);
            free(seen);
            (void)fclose(stream);
            return CO2_FORCING_TABLE_UNSUPPORTED_LATITUDE;
        }
        index = ((size_t)(year - first_year) * CO2_FORCING_MONTHS +
                 (size_t)(month - 1)) * CO2_FORCING_LATITUDE_BANDS +
                (size_t)band;
        if (seen[index] != 0)
        {
            if (error != NULL)
            {
                error->year = year;
                error->month = month;
                error->band = band;
            }
            free(values);
            free(seen);
            (void)fclose(stream);
            return CO2_FORCING_TABLE_DUPLICATE_RECORD;
        }
        values[index] = value;
        seen[index] = 1;
    }
    if (fclose(stream) != 0)
    {
        free(values);
        free(seen);
        return CO2_FORCING_TABLE_CLOSE_FAILED;
    }
    for (size_t index = 0; index < value_count; ++index)
    {
        if (seen[index] == 0)
        {
            if (error != NULL)
            {
                error->first_year = first_year;
                error->last_year = last_year;
            }
            free(values);
            free(seen);
            return CO2_FORCING_TABLE_INCOMPLETE;
        }
    }
    free(seen);
    table->values = values;
    table->first_year = first_year;
    table->last_year = last_year;
    return CO2_FORCING_TABLE_OK;
}

const float *co2_forcing_table_lookup(
    const co2_forcing_table *table, int year, int month)
{
    size_t index;

    if (table == NULL || table->values == NULL ||
        month < 1 || month > CO2_FORCING_MONTHS ||
        year < table->first_year || year > table->last_year)
    {
        return NULL;
    }
    index = ((size_t)(year - table->first_year) * CO2_FORCING_MONTHS +
             (size_t)(month - 1)) * CO2_FORCING_LATITUDE_BANDS;
    return table->values + index;
}

void co2_forcing_table_free(co2_forcing_table *table)
{
    if (table != NULL)
    {
        free(table->values);
        table->values = NULL;
        table->first_year = 0;
        table->last_year = 0;
    }
}
