#ifndef PLASIC_SHT_KERNELS_H
#define PLASIC_SHT_KERNELS_H

#ifdef __cplusplus
extern "C" {
#endif

enum plasic_sht_status
{
    PLASIC_SHT_OK = 0,
    PLASIC_SHT_INVALID_ARGUMENT = 1,
    PLASIC_SHT_INITIALIZATION_FAILED = 2,
    PLASIC_SHT_ALLOCATION_FAILED = 3,
    PLASIC_SHT_MPI_FAILED = 4
};

/*
 * Initialize the SHTns-backed spherical transform plan. Each MPI rank owns
 * nlpp consecutive Gaussian latitudes beginning at latitude_offset.
 */
int sht_init(int ntru, int nlat, int nlon, int nlpp, int latitude_offset);
void sht_finalize(void);

/* Scalar synthesis and analysis. Spectra use PlaSiC's interleaved
 * real/imaginary triangular layout; grid fields use local latitude rows. */
int sht_scalar_to_grid(int levels, const float *spectral, float *grid);
/* Gather a rank-local spectrum, then synthesize this rank's grid slice.
 * local_spectral_values must be ceil(global_values / mp_size()), i.e. the
 * per-rank block size; the last rank may hold fewer valid coefficients, and
 * the gather communicates only those. */
int sht_local_scalar_to_grid(
    int levels, const float *local_spectral, int local_spectral_values,
    float *grid);
int sht_grid_to_scalar(int levels, const float *grid, float *spectral);

/* Robert-form gradient of a scalar: both outputs are multiplied by cos(lat).
 * The first output points northward and the second eastward. */
int sht_scalar_to_grid_gradient(
    const float *spectral, float *northward, float *eastward);

/* Vector synthesis/analysis between divergence-vorticity spectra and
 * Robert-form zonal/meridional wind components. */
int sht_vortdiv_to_wind(
    int levels, float planetary_vorticity, const float *divergence,
    const float *vorticity, float *zonal_wind, float *meridional_wind);
int sht_wind_to_vortdiv(
    int levels, const float *zonal_wind, const float *meridional_wind,
    float *divergence, float *vorticity);

/* Adiabatic tendency projections. Inputs are local grid-point fields;
 * outputs are complete PlaSiC spectra on every rank. */
int sht_qtend(
    int levels, const float *scalar_source, const float *zonal_flux,
    const float *meridional_flux, float *tendency);
int sht_mktend(
    int levels, const float *temperature_source,
    const float *zonal_vorticity_flux,
    const float *meridional_vorticity_flux,
    const float *kinetic_plus_geopotential,
    const float *zonal_temperature_flux,
    const float *meridional_temperature_flux, float *divergence_tendency,
    float *temperature_tendency, float *vorticity_tendency);

#ifdef __cplusplus
}
#endif

#endif
