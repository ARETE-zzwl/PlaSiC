#include "calmod_kernels.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

/*
 * CALMOD 使用明确的 C 类型 int32_t / int64_t 和 float。
 *
 * 全模式只有本文件实现"步数 ↔ 日期"的换算，物理过程、输出和月平均
 * NetCDF 都通过这里的接口读取日历信息，不再各自复制日历规则。
 *
 * Calendar arithmetic lives only in this file: the physics, output and
 * monthly-NetCDF modules call these interfaces instead of re-implementing
 * calendar rules.
 */

enum
{
    MONTHS_PER_YEAR = 12,
    MINUTES_PER_DAY = 1440,
    MINUTES_PER_HOUR = 60,
    HOURS_PER_DAY = 24,
    /* 第 0 天的星期编号；仅用于兼容历史输出，不参与物理。 */
    WEEKDAY_OFFSET = 5
};

/*
 * Gregorian 平年每月天数，下标 0 为占位元素，使下标可直接使用月份编号。
 * 闰年的 2 月由 days_in_month() 单独加一天。
 */
static const int32_t gregorian_month_lengths[MONTHS_PER_YEAR + 1] = {
    0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

/* 每个模型日包含的时间步数，由 calmod_init() 写入。 */
static int32_t calendar_timesteps_per_day = 0;

static int32_t gregorian_is_leap_year(int32_t year)
{
    return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

/*
 * 第 0 年 1 月 1 日到 year 年 1 月 1 日的天数（前推 Gregorian）。
 * 公式按 400/100/4 年周期累加闰日，year 必须非负。
 */
static int64_t gregorian_days_before_year(int32_t year)
{
    const int64_t y = (int64_t)year;

    return 365 * y + (y + 3) / 4 - (y + 99) / 100 + (y + 399) / 400;
}

/*
 * 由"从纪元起算的完整天数"反解 Gregorian 年份与年内日偏移。
 *
 * 先用 400 年周期的平均年长估算年份，再用 days_before_year() 修正；
 * 估算误差不超过一年，因此循环最多执行两三次。
 */
static void gregorian_date_from_days(
    int64_t all_days, int32_t *year, int32_t *yday0)
{
    /* 400 年周期 = 146097 天，用平均年长估算年份，误差不超过一年。 */
    int32_t candidate = (int32_t)((all_days * 400) / 146097);

    while (gregorian_days_before_year(candidate) > all_days)
    {
        --candidate;
    }
    while (gregorian_days_before_year(candidate + 1) <= all_days)
    {
        ++candidate;
    }

    *year = candidate;
    *yday0 = (int32_t)(all_days - gregorian_days_before_year(candidate));
}

/* 由年内日偏移（从 0 开始）得到月份和月内日期（从 1 开始）。 */
static void gregorian_month_day(
    int32_t year, int32_t yday0, int32_t *month, int32_t *day)
{
    int32_t remaining = yday0;
    int32_t current = 1;

    while (remaining >=
           days_in_month(PLASIC_CALENDAR_GREGORIAN, year, current))
    {
        remaining -= days_in_month(PLASIC_CALENDAR_GREGORIAN, year, current);
        ++current;
    }

    *month = current;
    *day = remaining + 1;
}

int32_t calmod_init(int32_t timesteps_per_day)
{
    if (timesteps_per_day <= 0)
    {
        return CALMOD_INVALID_DATE;
    }
    calendar_timesteps_per_day = timesteps_per_day;
    return CALMOD_OK;
}

int32_t calendar_is_valid(int32_t calendar)
{
    return calendar == PLASIC_CALENDAR_GREGORIAN ||
           calendar == PLASIC_CALENDAR_FIXED_360;
}

const char *calmod_calendar_name(int32_t calendar)
{
    switch (calendar)
    {
    case PLASIC_CALENDAR_GREGORIAN:
        return "gregorian";
    case PLASIC_CALENDAR_FIXED_360:
        return "360day";
    default:
        return "unknown";
    }
}

int32_t calmod_calendar_from_name(const char *name)
{
    if (name == NULL)
    {
        return -1;
    }
    if (strcmp(name, "gregorian") == 0)
    {
        return PLASIC_CALENDAR_GREGORIAN;
    }
    if (strcmp(name, "360day") == 0 || strcmp(name, "360_day") == 0)
    {
        return PLASIC_CALENDAR_FIXED_360;
    }
    return -1;
}

const char *calmod_cf_calendar_name(int32_t calendar)
{
    return calendar == PLASIC_CALENDAR_FIXED_360 ? "360_day"
                                                 : "proleptic_gregorian";
}

const char *calmod_status_text(int32_t status)
{
    switch (status)
    {
    case CALMOD_OK:
        return "ok";
    case CALMOD_INVALID_CALENDAR:
        return "unknown calendar";
    case CALMOD_INVALID_DATE:
        return "invalid date or time";
    case CALMOD_TIME_NOT_ON_STEP:
        return "time does not fall on a time-step boundary";
    case CALMOD_STEP_OUT_OF_RANGE:
        return "date exceeds the representable step range";
    case CALMOD_NOT_INITIALIZED:
        return "calendar was not initialized";
    default:
        return "unknown calendar error";
    }
}

int32_t days_in_year(int32_t calendar, int32_t year)
{
    if (calendar == PLASIC_CALENDAR_FIXED_360)
    {
        return 360;
    }
    return gregorian_is_leap_year(year) ? 366 : 365;
}

int32_t days_in_month(int32_t calendar, int32_t year, int32_t month)
{
    if (calendar == PLASIC_CALENDAR_FIXED_360)
    {
        return month >= 1 && month <= MONTHS_PER_YEAR ? 30 : 0;
    }
    if (month < 1 || month > MONTHS_PER_YEAR)
    {
        return 0;
    }
    if (month == 2 && gregorian_is_leap_year(year))
    {
        return 29;
    }
    return gregorian_month_lengths[month];
}

int64_t days_before_year(int32_t calendar, int32_t year)
{
    if (calendar == PLASIC_CALENDAR_FIXED_360)
    {
        return 360LL * (int64_t)year;
    }
    return gregorian_days_before_year(year);
}

int64_t days_before_month(int32_t calendar, int32_t year, int32_t month)
{
    int64_t days;
    int32_t current;

    if (month < 1)
    {
        return days_before_year(calendar, year);
    }
    if (month > MONTHS_PER_YEAR)
    {
        return days_before_year(calendar, year) +
               (int64_t)days_in_year(calendar, year);
    }
    if (calendar == PLASIC_CALENDAR_FIXED_360)
    {
        return 360LL * (int64_t)year + 30LL * (int64_t)(month - 1);
    }

    days = gregorian_days_before_year(year);
    for (current = 1; current < month; ++current)
    {
        days += days_in_month(calendar, year, current);
    }
    return days;
}

int32_t cal2step(
    int32_t calendar, int32_t year, int32_t yday, int32_t hour,
    int32_t minute, int32_t *step)
{
    int64_t days;
    int64_t minute_of_day;
    int64_t total;

    if (!calendar_is_valid(calendar))
    {
        return CALMOD_INVALID_CALENDAR;
    }
    if (step == NULL || year < 0 || hour < 0 || hour >= HOURS_PER_DAY ||
        minute < 0 || minute >= MINUTES_PER_HOUR)
    {
        return CALMOD_INVALID_DATE;
    }
    if (yday < 1 || yday > days_in_year(calendar, year))
    {
        return CALMOD_INVALID_DATE;
    }
    if (calendar_timesteps_per_day <= 0)
    {
        return CALMOD_NOT_INITIALIZED;
    }

    days = days_before_year(calendar, year) + (int64_t)(yday - 1);
    minute_of_day = (int64_t)hour * MINUTES_PER_HOUR + (int64_t)minute;

    /*
     * 日内时刻必须精确落在某个时间步边界上：
     * 一天有 ntspd 个步，一步对应 1440/ntspd 分钟。
     */
    if ((minute_of_day * calendar_timesteps_per_day) % MINUTES_PER_DAY != 0)
    {
        return CALMOD_TIME_NOT_ON_STEP;
    }

    total = days * calendar_timesteps_per_day +
            (minute_of_day * calendar_timesteps_per_day) / MINUTES_PER_DAY;
    if (total < 0 || total > INT32_MAX)
    {
        return CALMOD_STEP_OUT_OF_RANGE;
    }

    *step = (int32_t)total;
    return CALMOD_OK;
}

void step2cal(int32_t calendar, int32_t step, int32_t datim[7])
{
    int64_t all_days;
    int32_t step_in_day;
    int32_t year;
    int32_t month;
    int32_t day;
    int32_t yday0;
    int32_t minute_of_day;

    if (datim == NULL)
    {
        return;
    }
    if (!calendar_is_valid(calendar) || step < 0 ||
        calendar_timesteps_per_day <= 0)
    {
        for (year = 0; year < 7; ++year)
        {
            datim[year] = -1;
        }
        return;
    }

    all_days = (int64_t)step / calendar_timesteps_per_day;
    step_in_day = (int32_t)((int64_t)step % calendar_timesteps_per_day);

    if (calendar == PLASIC_CALENDAR_FIXED_360)
    {
        year = (int32_t)(all_days / 360);
        yday0 = (int32_t)(all_days % 360);
        month = yday0 / 30 + 1;
        day = yday0 % 30 + 1;
    }
    else
    {
        gregorian_date_from_days(all_days, &year, &yday0);
        gregorian_month_day(year, yday0, &month, &day);
    }

    minute_of_day = (int32_t)(
        ((int64_t)step_in_day * MINUTES_PER_DAY) /
        calendar_timesteps_per_day);

    datim[0] = year;
    datim[1] = month;
    datim[2] = day;
    datim[3] = minute_of_day / MINUTES_PER_HOUR;
    datim[4] = minute_of_day % MINUTES_PER_HOUR;
    datim[5] = (int32_t)((all_days + WEEKDAY_OFFSET) % 7);
    datim[6] = calendar == PLASIC_CALENDAR_GREGORIAN &&
                       gregorian_is_leap_year(year)
                   ? 1
                   : 0;
}

int32_t ndayofyear(int32_t calendar, int32_t step)
{
    int32_t datim[7];
    int64_t days_before_month_in_year;

    step2cal(calendar, step, datim);
    if (datim[0] < 0)
    {
        return 0;
    }

    days_before_month_in_year =
        days_before_month(calendar, datim[0], datim[1]) -
        days_before_year(calendar, datim[0]);
    return datim[2] + (int32_t)days_before_month_in_year;
}

void momint(
    int32_t calendar, int32_t step, int32_t *month_a, int32_t *month_b,
    float *weight)
{
    int32_t datim[7] = {0, 0, 0, 0, 0, 0, 0};
    int32_t neighbor_month;
    int32_t neighbor_year;
    float fractional_day;
    float days_month_a;
    float days_month_b;
    float middle_month_a;
    float middle_month_b;

    if (month_a == NULL || month_b == NULL || weight == NULL)
    {
        return;
    }

    step2cal(calendar, step, datim);
    if (datim[0] < 0)
    {
        *month_a = -1;
        *month_b = -1;
        *weight = 0.0f;
        return;
    }

    *month_a = datim[1];
    fractional_day = (float)datim[2] +
                     (((float)datim[3] * 60.0f + (float)datim[4]) /
                      (float)MINUTES_PER_DAY);

    days_month_a = (float)days_in_month(calendar, datim[0], *month_a);
    middle_month_a = 0.5f * (days_month_a + 1.0f);

    if (fractional_day > middle_month_a)
    {
        *month_b = *month_a + 1;
        neighbor_month = *month_b;
        neighbor_year = datim[0];
        if (neighbor_month > MONTHS_PER_YEAR)
        {
            neighbor_month = 1;
            ++neighbor_year;
        }
    }
    else
    {
        *month_b = *month_a - 1;
        neighbor_month = *month_b;
        neighbor_year = datim[0];
        if (neighbor_month < 1)
        {
            neighbor_month = MONTHS_PER_YEAR;
            --neighbor_year;
        }
    }

    /*
     * 相邻月份只会落在 12 月或 1 月，月长与年份无关；把年份限制在
     * 非负范围只是为了让 days_in_month() 的前置条件始终成立。
     */
    if (neighbor_year < 0)
    {
        neighbor_year = 0;
    }
    days_month_b = (float)days_in_month(calendar, neighbor_year, neighbor_month);
    middle_month_b = 0.5f * (days_month_b + 1.0f);

    *weight = fabsf(fractional_day - middle_month_a) /
              (middle_month_a + middle_month_b - 1.0f);
}
