#ifndef PLASIC_MONTHLY_NETCDF_H
#define PLASIC_MONTHLY_NETCDF_H

int runtime_monthly_netcdf_open(const char *directory);
int runtime_monthly_netcdf_accumulate(int year, int month);
int runtime_monthly_netcdf_close(void);
void runtime_monthly_netcdf_abort(void);

#endif
