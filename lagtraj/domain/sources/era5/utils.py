"""
ECMWF ERA5 specific functionality needed for trajectory integration
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import datetime

from .... import njit
from .constants import rg, rd, rv_over_rd_minus_one


def _load_ecmwf_level_coefficients():
    levels_file = Path(__file__).parent / "137levels.dat"
    levels_table = pd.read_table(levels_file, sep=r"\s+")

    ds_levels_coeffs = xr.Dataset(coords=dict(level=levels_table["n"].values),)
    ds_levels_coeffs["a"] = xr.DataArray(levels_table["a[Pa]"].values, dims=("level",))
    ds_levels_coeffs["b"] = xr.DataArray(levels_table["b"].values, dims=("level",))
    return ds_levels_coeffs


ds_levels_coeffs = _load_ecmwf_level_coefficients()


@njit
def _calculate_heights_and_pressures(
    p_surf, height_surf, a_coeffs, b_coeffs, t_field, q_field
):
    """Calculate heights and pressures at model levels using
    the hydrostatic equation (not taking into account hydrometeors).
    """
    k_max = t_field.shape[0]
    j_max = t_field.shape[1]
    i_max = t_field.shape[2]
    height_h = np.empty((k_max, j_max, i_max))
    height_f = np.empty((k_max, j_max, i_max))
    p_h = np.empty((k_max, j_max, i_max))
    p_f = np.empty((k_max, j_max, i_max))
    rd_over_rg = rd / rg
    for i in range(i_max):
        for j in range(j_max):
            p_s = p_surf[j, i]
            p_h[k_max - 1, j, i] = p_s
            p_h_k_plus = p_s
            z_s = height_surf[j, i]
            height_h[k_max - 1, j, i] = z_s
            height_h_k_plus = z_s
            for k in range(k_max - 2, -1, -1):
                # Pressure at this half level
                p_h_k = a_coeffs[k] + (b_coeffs[k] * p_s)
                p_h[k, j, i] = p_h_k
                # Pressure at corresponding full level
                p_f_k_plus = 0.5 * (p_h_k + p_h_k_plus)
                p_f[k + 1, j, i] = p_f_k_plus
                # Virtual temperature
                tvrd_over_rg = (
                    rd_over_rg
                    * t_field[k, j, i]
                    * (1.0 + rv_over_rd_minus_one * q_field[k, j, i])
                )
                # Integration to half level first
                height_f[k + 1, j, i] = height_h_k_plus + (
                    tvrd_over_rg * np.log(p_h_k_plus / p_f_k_plus)
                )
                # Integration to full levels
                # reset of scalar temporary variables
                height_h_k_plus = height_h_k_plus + (
                    tvrd_over_rg * np.log(p_h_k_plus / p_h_k)
                )
                height_h[k, j, i] = height_h_k_plus
                p_h_k_plus = p_h_k
            p_f_k_plus = 0.5 * p_h_k_plus
            p_f[0, j, i] = p_f_k_plus
            height_f[0, j, i] = height_h_k_plus + (
                tvrd_over_rg * np.log(p_h_k_plus / p_f_k_plus)
            )

    return height_h, height_f, p_h, p_f


def calculate_heights_and_pressures(ds):
    """Calculates height and pressure fields to ERA5 model level data arrays"""

    datasets = []
    ds_ = ds
    # `calculate_heights_and_pressures` expends the data to have the shape
    # (level, lat, lon) so we explicitly expand the dataset here and squeeze
    # later
    required_dims = ("time", "level", "lat", "lon")
    missing_dims = list(filter(lambda d: d not in ds.dims, required_dims))
    ds_ = ds_.expand_dims(missing_dims).transpose(*required_dims)

    levels = ds.level
    if levels.values[0] > levels.values[-1]:
        raise Exception(
            "Height and pressure calculation assumes top-down " "ordering of the levels"
        )

    a_coeffs = ds_levels_coeffs.sel(level=ds.level).a.values
    b_coeffs = ds_levels_coeffs.sel(level=ds.level).b.values

    for t in ds_.time.values:
        ds_time = ds_.sel(time=t)
        p_surf = ds_time.sp.values
        # Convert from geopotential to height
        height_surf = ds_time.z.values / rg
        t_field = ds_time.t.values
        q_field = ds_time.q.values

        height_dims = ds_time.t.dims
        height_h, height_f, p_h, p_f = _calculate_heights_and_pressures(
            p_surf, height_surf, a_coeffs, b_coeffs, t_field, q_field,
        )
        ds_extra = xr.Dataset(coords=ds_time.coords)
        ds_extra["height_h"] = xr.DataArray(
            height_h,
            dims=height_dims,
            attrs={"long_name": "height above sea level at half level", "units": "m",},
        )
        ds_extra["height_f"] = xr.DataArray(
            height_f,
            dims=height_dims,
            attrs={"long_name": "height above sea level at full level", "units": "m",},
        )
        ds_extra["p_h"] = xr.DataArray(
            p_h,
            dims=height_dims,
            attrs={"long_name": "pressure at half level", "units": "Pa"},
        )
        ds_extra["p_f"] = xr.DataArray(
            p_f,
            dims=height_dims,
            attrs={"long_name": "pressure at full level", "units": "Pa"},
        )
        datasets.append(ds_extra)

    return xr.concat(datasets, dim="time").squeeze()


def add_era5_global_attributes(ds, creation_datetime):
    """Adds global attributes to datasets"""
    global_attrs = {
        r"conventions": r"CF-1.7",
        r"contact": r"l.c.denby[at]leeds[dot]ac[dot again]uk s.boeing[at]leeds[dot]ac[dot again]uk",
        r"era5_reference": r"Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., ... & Simmons, A. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society.",
        r"created": creation_datetime.isoformat(),
        r"created_with": r"https://github.com/EUREC4A-UK/lagtraj",
        r"note": "Contains modified Copernicus Service information ",
    }
    for attribute in global_attrs:
        ds.attrs[attribute] = global_attrs[attribute]
