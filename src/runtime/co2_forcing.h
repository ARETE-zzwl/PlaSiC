#ifndef PLASIC_CO2_FORCING_H
#define PLASIC_CO2_FORCING_H

int runtime_co2_forcing_open(const char *path, float fixed_ppm);
int runtime_co2_forcing_apply(int year, int month);
void runtime_co2_forcing_close(void);
const char *runtime_co2_forcing_source(void);

#endif
