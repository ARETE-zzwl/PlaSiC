#ifndef PLASIC_CO2_FORCING_TABLE_H
#define PLASIC_CO2_FORCING_TABLE_H

#include <stdint.h>

enum co2_forcing_table_status
{
    CO2_FORCING_TABLE_OK = 0,
    CO2_FORCING_TABLE_OPEN_FAILED = 1,
    CO2_FORCING_TABLE_BAD_HEADER = 2,
    CO2_FORCING_TABLE_BAD_RECORD = 3,
    CO2_FORCING_TABLE_UNSUPPORTED_LATITUDE = 4,
    CO2_FORCING_TABLE_DUPLICATE_RECORD = 5,
    CO2_FORCING_TABLE_INCOMPLETE = 6,
    CO2_FORCING_TABLE_TOO_LARGE = 7,
    CO2_FORCING_TABLE_ALLOCATION_FAILED = 8,
    CO2_FORCING_TABLE_SCAN_FAILED = 9,
    CO2_FORCING_TABLE_REWIND_FAILED = 10,
    CO2_FORCING_TABLE_CLOSE_FAILED = 11
};

enum
{
    CO2_FORCING_MONTHS = 12,
    CO2_FORCING_LATITUDE_BANDS = 12
};

/* Machine-readable details for the status codes above. */
typedef struct co2_forcing_table_error
{
    int32_t year;
    int32_t month;
    int32_t band;
    float latitude;
    int32_t first_year;
    int32_t last_year;
} co2_forcing_table_error;

/* Complete serial representation of a CO2 forcing file. */
typedef struct co2_forcing_table
{
    float *values;
    int32_t first_year;
    int32_t last_year;
} co2_forcing_table;

/*
 * Read a CMIP-style CO2 forcing CSV and retain its complete
 * [year][month][latitude band] table. The reader performs no parallel
 * communication and no runtime error reporting: on failure the table stays
 * empty and error carries the details needed to build a diagnostic message.
 * On success the table owns values until co2_forcing_table_free() is called.
 */
int co2_forcing_table_load(
    const char *path, co2_forcing_table *table,
    co2_forcing_table_error *error);

/* Return the 12 latitude bands for one month, or NULL when unavailable. */
const float *co2_forcing_table_lookup(
    const co2_forcing_table *table, int year, int month);

void co2_forcing_table_free(co2_forcing_table *table);

#endif
