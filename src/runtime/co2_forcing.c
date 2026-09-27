/*
 * CO2 强迫的运行时封装：把命令行选择的 CO2 模式装载到物理模块。
 * Runtime wrapper for CO2 forcing: loads the CO2 mode selected on the command
 * line into the physics module.
 *
 * 两种互斥模式 / Two mutually exclusive modes:
 *   1. 固定浓度：--co2-ppm PPMV，所有纬度带使用同一浓度值。
 *      Fixed concentration: --co2-ppm PPMV applies one value to all bands.
 *   2. 数据表：--co2-forcing FILE，按月读取 CMIP 风格 CSV，索引为
 *      [year][month][latitude band]。
 *      Table: --co2-forcing FILE reads a CMIP-style CSV indexed by
 *      [year][month][latitude band].
 */

#include "runtime_internal.h"

#include "co2_forcing.h"
#include "tools/co2_forcing_table.h"

typedef struct co2_forcing_state
{
    /* 完整的 [year][month][band] 数据表；固定浓度模式下为空。 */
    /* Complete [year][month][band] table; empty in fixed-concentration mode.
     */
    co2_forcing_table table;

    /* 面向诊断输出的来源描述：CSV 路径或固定浓度文本。 */
    /* Source description for diagnostics: CSV path or fixed-concentration text. */
    char source[4096];

    /* 固定浓度值（ppmv）；<= 0 表示当前处于数据表模式。 */
    /* Fixed concentration in ppmv; <= 0 means table mode is active. */
    float fixed_ppm;

    /* 最近一次成功 apply 的模型年月，用于跳过重复映射。 */
    /* Model date of the last successful apply, used to skip redundant mapping. */
    int32_t active_year;
    int32_t active_month;
} co2_forcing_state;

static co2_forcing_state forcing;

/*
 * Map a table-layer status code to a runtime status code.
 */
static int co2_table_status_to_runtime(int status)
{
    return status == CO2_FORCING_TABLE_ALLOCATION_FAILED
               ? PLASIC_RUNTIME_INITIALIZATION_FAILED
               : PLASIC_RUNTIME_IO_FAILED;
}

/*
 * Translate a failed table load into a contextual error message on the root
 * rank.
 * Only the root rank calls this because only it ran co2_forcing_table_load().
 */
static void report_co2_table_error(
    const char *path, int status, const co2_forcing_table_error *error)
{
    switch (status)
    {
    case CO2_FORCING_TABLE_OPEN_FAILED:
        runtime_set_error("cannot open CO2 forcing file '%s'", path);
        break;
    case CO2_FORCING_TABLE_BAD_HEADER:
        runtime_set_error("invalid CO2 forcing CSV header in '%s'", path);
        break;
    case CO2_FORCING_TABLE_BAD_RECORD:
        runtime_set_error("invalid CO2 forcing record in '%s'", path);
        break;
    case CO2_FORCING_TABLE_UNSUPPORTED_LATITUDE:
        runtime_set_error(
            "unsupported CO2 forcing latitude %.6g", (double)error->latitude);
        break;
    case CO2_FORCING_TABLE_DUPLICATE_RECORD:
        runtime_set_error(
            "duplicate CO2 forcing record for %d-%02d band %d",
            error->year, error->month, error->band);
        break;
    case CO2_FORCING_TABLE_INCOMPLETE:
        runtime_set_error(
            "CO2 forcing file is incomplete for years %d-%d",
            error->first_year, error->last_year);
        break;
    case CO2_FORCING_TABLE_TOO_LARGE:
        runtime_set_error("CO2 forcing table is too large");
        break;
    case CO2_FORCING_TABLE_ALLOCATION_FAILED:
        runtime_set_error("could not allocate CO2 forcing table");
        break;
    case CO2_FORCING_TABLE_SCAN_FAILED:
        runtime_set_error("could not scan CO2 forcing file '%s'", path);
        break;
    case CO2_FORCING_TABLE_REWIND_FAILED:
        runtime_set_error("could not rewind CO2 forcing file '%s'", path);
        break;
    case CO2_FORCING_TABLE_CLOSE_FAILED:
        runtime_set_error("could not close CO2 forcing file '%s'", path);
        break;
    default:
        runtime_set_error("could not read CO2 forcing file '%s'", path);
        break;
    }
}

/*
 * Open CO2 forcing
 *
 * path      非空 CSV 路径进入数据表模式；否则必须提供 fixed_ppm。
 *           A non-empty CSV path selects table mode; otherwise fixed_ppm is
 *           required.
 */
int runtime_co2_forcing_open(const char *path, float fixed_ppm)
{
    /* metadata 依次为：加载状态、起始年、结束年。 */
    /* metadata holds: load status, first year, last year, in that order. */
    int32_t metadata[3] = {PLASIC_RUNTIME_OK, 0, 0};

    co2_forcing_table_error error;

    /* 用于 MPI 归约的分配失败标志：0 表示各进程都分配成功。 */
    /* Allocation-failure flag reduced over MPI: 0 means every rank succeeded.
     */
    float allocation_failure = 0.0f;

    size_t value_count;
    int table_status = CO2_FORCING_TABLE_OK;

    const int has_path = path != NULL && path[0] != '\0';
    const int has_fixed = isfinite(fixed_ppm) && fixed_ppm > 0.0f;

    /* 先释放上次打开的强迫 */
    /* Release any previous forcing first */
    runtime_co2_forcing_close();

    if (!isfinite(fixed_ppm) || fixed_ppm < 0.0f ||
        has_path + has_fixed != 1)
    {
        runtime_set_error(
            "specify exactly one CO2 mode: --co2-forcing FILE or "
            "--co2-ppm PPMV");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    if (has_fixed)
    {
        forcing.fixed_ppm = fixed_ppm;
        (void)snprintf(
            forcing.source, sizeof(forcing.source),
            "fixed CO2: %.9g ppmv", (double)fixed_ppm);
        forcing.active_year = -1;
        forcing.active_month = -1;
        return PLASIC_RUNTIME_OK;
    }

    if (strlen(path) >= sizeof(forcing.source))
    {
        runtime_set_error("a valid CO2 forcing CSV path is required");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    if (mp_is_root())
    {
        /* 只让根进程做一次文件扫描与解析，随后广播结果。 */
        /* Only the root rank scans and parses the file; results are broadcast afterwards. */
        table_status = co2_forcing_table_load(path, &forcing.table, &error);
        metadata[0] = table_status == CO2_FORCING_TABLE_OK
                          ? PLASIC_RUNTIME_OK
                          : co2_table_status_to_runtime(table_status);
        metadata[1] = forcing.table.first_year;
        metadata[2] = forcing.table.last_year;
    }

    /* 先把状态和时间范围同步给所有进程，让它们可以决定后续行为。 */
    /* Synchronize status and time range first so every rank can decide what to do next. */
    if (mp_broadcast_integer(metadata, 3) != 0)
    {
        runtime_set_error("could not broadcast CO2 forcing metadata");
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }

    if (metadata[0] != PLASIC_RUNTIME_OK)
    {
        if (mp_is_root())
        {
            report_co2_table_error(path, table_status, &error);
        }
        return (int)metadata[0];
    }

    /* 所有进程统一重建首末年份，并据此算出完整表长度。 */
    /* Rebuild the year range consistently on every rank and derive the full table length from it. */
    forcing.table.first_year = metadata[1];
    forcing.table.last_year = metadata[2];
    value_count =
        ((size_t)(forcing.table.last_year - forcing.table.first_year) + 1) *
        CO2_FORCING_MONTHS * CO2_FORCING_LATITUDE_BANDS;

    if (!mp_is_root())
    {
        /* 根进程的表由 co2_forcing_table_load() 分配，非根进程在此自行分配。
         */
        /* The root table was allocated by co2_forcing_table_load(); non-root ranks allocate theirs here. */
        forcing.table.values =
            (float *)malloc(value_count * sizeof(*forcing.table.values));
    }

    /* 用一次全局归约确认没有进程分配失败，避免部分进程继续运行。 */
    /* One global reduction confirms no rank failed to allocate, preventing split execution. */
    allocation_failure = forcing.table.values == NULL ? 1.0f : 0.0f;
    if (mp_allreduce_real(&allocation_failure, 1) != 0 ||
        allocation_failure != 0.0f)
    {
        runtime_set_error("could not allocate distributed CO2 forcing table");
        runtime_co2_forcing_close();
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }

    /* 把根进程读到的全部浓度值广播到每个进程的表缓冲区。 */
    /* Broadcast every concentration read by the root to each rank's table buffer. */
    if (mp_broadcast_real(forcing.table.values, (int32_t)value_count) != 0)
    {
        runtime_set_error("could not broadcast CO2 forcing values");
        runtime_co2_forcing_close();
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    (void)snprintf(forcing.source, sizeof(forcing.source), "%s", path);

    /* 复位缓存，确保第一次 apply 一定会真正映射到物理模块。 */
    /* Reset the cache so the first apply always maps to the physics module. */
    forcing.active_year = -1;
    forcing.active_month = -1;
    return PLASIC_RUNTIME_OK;
}

/*
 * 把给定模型年月的 CO2 纬度带浓度应用（映射）到物理模块。
 * Apply (map) the CO2 latitude-band concentrations of the given model date to
 * the physics module.
 */
int runtime_co2_forcing_apply(int year, int month)
{
    const float *bands;
    float fixed_bands[CO2_FORCING_LATITUDE_BANDS];

    if (forcing.fixed_ppm > 0.0f)
    {
        /* 固定浓度模式：所有纬度带使用相同值。 */
        /* Fixed-concentration mode: every latitude band gets the same value.
         */
        if (forcing.active_year == year && forcing.active_month == month)
        {
            return PLASIC_RUNTIME_OK;
        }
        for (int band = 0; band < CO2_FORCING_LATITUDE_BANDS; ++band)
        {
            fixed_bands[band] = forcing.fixed_ppm;
        }
        if (plasic_physics_set_co2_bands(
                fixed_bands, CO2_FORCING_LATITUDE_BANDS) != PLASIC_RUNTIME_OK)
        {
            runtime_set_error("could not apply fixed CO2 concentration");
            return PLASIC_RUNTIME_INITIALIZATION_FAILED;
        }
        forcing.active_year = year;
        forcing.active_month = month;
        return PLASIC_RUNTIME_OK;
    }

    /* 数据表模式：先查表，日期超出表范围会得到 NULL。 */
    /* Table mode: look up first; dates outside the table range yield NULL. */
    bands = co2_forcing_table_lookup(&forcing.table, year, month);
    if (year == forcing.active_year && month == forcing.active_month)
    {
        return PLASIC_RUNTIME_OK;
    }
    if (bands == NULL)
    {
        runtime_set_error(
            "CO2 forcing table does not cover model date %d-%02d",
            year, month);
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    if (plasic_physics_set_co2_bands(
            bands, CO2_FORCING_LATITUDE_BANDS) != PLASIC_RUNTIME_OK)
    {
        runtime_set_error("could not map CO2 forcing for %d-%02d", year, month);
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }
    forcing.active_year = year;
    forcing.active_month = month;
    return PLASIC_RUNTIME_OK;
}

void runtime_co2_forcing_close(void)
{
    co2_forcing_table_free(&forcing.table);
    memset(&forcing, 0, sizeof(forcing));
    forcing.active_year = -1;
    forcing.active_month = -1;
}

/*
 * 返回当前 CO2 来源的描述字符串（固定浓度文本或 CSV 路径）。
 * Return the description of the active CO2 source (fixed-concentration text
 * or CSV path).
 */
const char *runtime_co2_forcing_source(void)
{
    return forcing.source;
}
