from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class RegionBounds:
    west: float = -60.0
    east: float = 60.0
    south: float = -30.0
    north: float = 30.0

    @property
    def crosses_dateline(self) -> bool:
        return self.west > self.east

    def validate(self) -> None:
        if not all(np.isfinite((self.west, self.east, self.south, self.north))):
            raise ValueError("Region bounds must be finite")
        if not (-180.0 <= self.west <= 180.0 and -180.0 <= self.east <= 180.0):
            raise ValueError("Longitude bounds must be within -180° and 180°")
        if not (-90.0 <= self.south < self.north <= 90.0):
            raise ValueError("South must be smaller than north within -90°..90°")
        if self.west == self.east:
            raise ValueError("West and east must differ; use -180° and 180° for the globe")


@dataclass(frozen=True)
class RegionData:
    values: NDArray[np.floating]
    latitudes: NDArray[np.float64]
    longitudes: NDArray[np.float64]
    display_longitudes: NDArray[np.float64]
    crosses_dateline: bool


@lru_cache(maxsize=16)
def gaussian_grid(nlat: int, nlon: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Coordinates for the pstream Gaussian-latitude, regular-longitude grid."""
    if nlat < 2 or nlon < 2:
        raise ValueError("A plottable global grid needs at least 2×2 cells")
    nodes, _weights = np.polynomial.legendre.leggauss(nlat)
    latitudes = np.degrees(np.arcsin(nodes))[::-1]
    longitudes = np.arange(nlon, dtype=np.float64) * (360.0 / nlon)
    longitudes = (longitudes + 180.0) % 360.0 - 180.0
    return latitudes, longitudes


def region_subset(values: NDArray[np.floating], bounds: RegionBounds) -> RegionData:
    """Subset (..., latitude, longitude), preserving west-to-east order across the dateline."""
    bounds.validate()
    array = np.asanyarray(values)
    if array.ndim < 2:
        raise ValueError("Spatial data must include latitude and longitude dimensions")
    nlat, nlon = array.shape[-2:]
    latitudes, longitudes = gaussian_grid(nlat, nlon)

    lat_mask = (latitudes >= bounds.south) & (latitudes <= bounds.north)
    distance = (longitudes - bounds.west) % 360.0
    span = (bounds.east - bounds.west) % 360.0
    lon_mask = distance <= span + 1.0e-10
    lat_indices = np.flatnonzero(lat_mask)
    lon_indices = np.flatnonzero(lon_mask)
    if lat_indices.size < 2 or lon_indices.size < 2:
        raise ValueError("Selected region is empty or smaller than two grid cells")

    order = np.argsort(distance[lon_indices])
    lon_indices = lon_indices[order]
    display_longitudes = bounds.west + distance[lon_indices]
    subset = np.take(np.take(array, lat_indices, axis=-2), lon_indices, axis=-1)
    return RegionData(
        values=subset,
        latitudes=latitudes[lat_indices],
        longitudes=longitudes[lon_indices],
        display_longitudes=display_longitudes,
        crosses_dateline=bounds.crosses_dateline,
    )


def finite_color_range(
    values: Iterable[NDArray[np.floating]] | NDArray[np.floating],
    *,
    nonnegative: bool = False,
) -> tuple[float, float]:
    """Return a safe finite range, including deterministic all-missing/constant fallbacks."""
    arrays: Iterable[NDArray[np.floating]]
    if isinstance(values, np.ndarray) or np.ma.isMaskedArray(values):
        arrays = (values,)
    else:
        arrays = values
    lower = np.inf
    upper = -np.inf
    found = False
    for value in arrays:
        data = np.asanyarray(value)
        if np.ma.isMaskedArray(data):
            data = data.filled(np.nan)
        finite = np.asarray(data, dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        if finite.size:
            lower = min(lower, float(np.min(finite)))
            upper = max(upper, float(np.max(finite)))
            found = True
    if not found:
        return 0.0, 1.0
    if lower == upper:
        if nonnegative and lower == 0.0:
            return 0.0, 1.0
        padding = max(abs(lower) * 0.01, 0.5)
        padded_lower = lower - padding
        if nonnegative and lower >= 0.0:
            padded_lower = max(0.0, padded_lower)
        return padded_lower, upper + padding
    return lower, upper


def gaussian_global_mean(values: NDArray[np.floating]) -> tuple[float, float]:
    """Gaussian-quadrature global mean and finite-data area coverage fraction."""
    data = np.asanyarray(values)
    if data.ndim != 2:
        raise ValueError("Global mean requires a two-dimensional latitude-longitude field")
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    data = np.asarray(data, dtype=np.float64)
    _nodes, latitude_weights = np.polynomial.legendre.leggauss(data.shape[0])
    cell_weights = np.broadcast_to(latitude_weights[:, None], data.shape)
    valid = np.isfinite(data) & np.isfinite(cell_weights) & (cell_weights > 0.0)
    valid_weight = float(np.sum(cell_weights[valid], dtype=np.float64))
    total_weight = float(np.sum(cell_weights, dtype=np.float64))
    coverage = valid_weight / total_weight if total_weight > 0.0 else 0.0
    if valid_weight == 0.0:
        return float("nan"), coverage
    mean = float(np.sum(data[valid] * cell_weights[valid], dtype=np.float64) / valid_weight)
    return mean, coverage


def frame_time_label(frame: object) -> str:
    year = int(getattr(frame, "year"))
    month = int(getattr(frame, "month"))
    day = int(getattr(frame, "day"))
    hour = int(getattr(frame, "hour"))
    minute = int(getattr(frame, "minute"))
    step = int(getattr(frame, "step"))
    return f"{year:04d}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}  ·  step {step}"


def resample_latitude_for_image(values: NDArray[np.floating], latitudes: NDArray[np.float64], rows: int | None = None) -> NDArray[np.float32]:
    """Resample Gaussian rows to a linear-latitude display raster; statistics stay untouched."""
    data = np.asanyarray(values)
    if data.ndim != 2:
        raise ValueError("Display resampling requires a 2-D field")
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    output_rows = rows or data.shape[0]
    target = np.linspace(float(latitudes[-1]), float(latitudes[0]), output_rows)
    result = np.empty((output_rows, data.shape[1]), dtype=np.float32)
    source_x = latitudes[::-1]
    for column in range(data.shape[1]):
        source_y = np.asarray(data[::-1, column], dtype=np.float64)
        valid = np.isfinite(source_x) & np.isfinite(source_y)
        if np.count_nonzero(valid) < 2:
            result[:, column] = np.nan
        else:
            valid_x = source_x[valid]
            valid_y = source_y[valid]
            result[:, column] = np.interp(
                target,
                valid_x,
                valid_y,
                left=valid_y[0],
                right=valid_y[-1],
            )
    return result
