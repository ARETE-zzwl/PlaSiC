#ifndef PLASIC_CALMOD_KERNELS_H
#define PLASIC_CALMOD_KERNELS_H

#include <stdint.h>

/*
 * CALMOD：PlaSiC 唯一的日历算术实现。
 *
 * 日历类型决定 nstep 与日期的映射：
 *
 *   PLASIC_CALENDAR_GREGORIAN = 0
 *       前推 Gregorian 公历：平年 365 天、闰年 366 天，使用公历月长。
 *       第 0 年按 proleptic Gregorian 规则是闰年（可被 400 整除）。
 *
 *   PLASIC_CALENDAR_FIXED_360 = 1
 *       12 月 × 30 天的理想化日历：每年 360 天，无闰年。
 *
 * 两种日历共用同一纪元：
 *
 *       nstep = 0  ≡  第 0 年 1 月 1 日 00:00
 *
 * 日历类型写入 restart 文件，热启动时据此解释 nstep。
 */
enum plasic_calendar_kind
{
    PLASIC_CALENDAR_GREGORIAN = 0,
    PLASIC_CALENDAR_FIXED_360 = 1
};

/* 日历接口返回码。 */
enum calmod_status
{
    CALMOD_OK = 0,
    CALMOD_INVALID_CALENDAR = 1,
    CALMOD_INVALID_DATE = 2,
    CALMOD_TIME_NOT_ON_STEP = 3,
    CALMOD_STEP_OUT_OF_RANGE = 4,
    CALMOD_NOT_INITIALIZED = 5
};

/*
 * 写入每个模型日的时间步数，必须在调用任何换算接口之前执行一次。
 * timesteps_per_day <= 0 时返回 CALMOD_INVALID_DATE。
 */
int32_t calmod_init(int32_t timesteps_per_day);

/* 日历类型合法性检查；合法返回 1，否则返回 0。 */
int32_t calendar_is_valid(int32_t calendar);

/* CLI/JSON 名称与日历类型之间的转换。 */
const char *calmod_calendar_name(int32_t calendar);
int32_t calmod_calendar_from_name(const char *name);
/* CF/NetCDF 时间轴使用的日历名称。 */
const char *calmod_cf_calendar_name(int32_t calendar);

/* 返回码对应的说明文本。 */
const char *calmod_status_text(int32_t status);

/* 日期算术：年份必须非负，月份为 1–12。 */
int32_t days_in_year(int32_t calendar, int32_t year);
int32_t days_in_month(int32_t calendar, int32_t year, int32_t month);
int64_t days_before_year(int32_t calendar, int32_t year);
int64_t days_before_month(int32_t calendar, int32_t year, int32_t month);

/*
 * 日期时间 → 步数。
 *
 * year  年份，从 0 开始；
 * yday  年内日序，从 1 开始（Gregorian 下 ≤ 365/366，360 天下 ≤ 360）；
 * hour  时，0–23；
 * minute 分，0–59。
 *
 * 返回 CALMOD_OK 时写出 step；其余返回码表示参数错误，step 不被修改。
 * 若时/分不能精确落在某个时间步边界（即 (hour×60+minute)×ntspd 不能被
 * 1440 整除），返回 CALMOD_TIME_NOT_ON_STEP。
 */
int32_t cal2step(
    int32_t calendar, int32_t year, int32_t yday, int32_t hour,
    int32_t minute, int32_t *step);

/*
 * 步数 → 日期时间。
 *
 * datim[0..6] = 年、月、日、时、分、星期、闰年标志。
 * 360 天日历下星期按同一公式给出，闰年标志恒为 0。
 */
void step2cal(int32_t calendar, int32_t step, int32_t datim[7]);

/* 步数 → 年内日序（从 1 开始），供辐射季节循环使用。 */
int32_t ndayofyear(int32_t calendar, int32_t step);

/*
 * 步数 → 月平均气候场插值关系。
 *
 * month_a 为当前月份；month_b 为相邻月份（位于本月月中点之后取下一个月，
 * 否则取上一个月，因此可能为 0 或 13，按 14 槽气候场解释）；
 * weight 为当前时刻到 month_a 月中点的归一化距离。
 */
void momint(
    int32_t calendar, int32_t step, int32_t *month_a, int32_t *month_b,
    float *weight);

#endif
