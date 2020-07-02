"""
ECMWF ERA5 specific functionality needed for trajectory integration
"""
import os
import numbers

import numpy as np
import pandas as pd
import xarray as xr

from ..interpolation import steffen_3d, njit
from .constants import rg, rd, rv_over_rd_minus_one

levels_file = os.path.dirname(__file__) + "/137levels.dat"
levels_table = pd.read_table(levels_file, sep="\s+")
# skip the first row, as it corresponds to the top of the atmosphere
# which is not in the data
a_coeffs_137 = levels_table["a[Pa]"].values[1:]
b_coeffs_137 = levels_table["b"].values[1:]


@njit
def calculate_heights_and_pressures(
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


def add_heights_and_pressures(ds_from_era5):
    """Adds height and pressure fields to ERA5 model level data arrays"""
    len_temp = len(ds_from_era5["t"])
    shape_temp = np.shape(ds_from_era5["t"])
    ds_from_era5["height_h"] = (
        ("time", "level", "latitude", "longitude"),
        np.empty(shape_temp),
        {"long_name": "height above sea level at half level", "units": "metres"},
    )
    ds_from_era5["height_f"] = (
        ("time", "level", "latitude", "longitude"),
        np.empty(shape_temp),
        {"long_name": "height above sea level at full level", "units": "metres"},
    )
    ds_from_era5["p_h"] = (
        ("time", "level", "latitude", "longitude"),
        np.empty(shape_temp),
        {"long_name": "pressure at half level", "units": "Pa"},
    )
    ds_from_era5["p_f"] = (
        ("time", "level", "latitude", "longitude"),
        np.empty(shape_temp),
        {"long_name": "pressure at full level", "units": "Pa"},
    )
    for time_index in range(len_temp):
        p_surf = ds_from_era5["sp"].values[time_index, :, :]
        # Convert from geopotential to height
        height_surf = ds_from_era5["z"].values[time_index, :, :] / rg
        t_field = ds_from_era5["t"].values[time_index, :, :, :]
        q_field = ds_from_era5["q"].values[time_index, :, :, :]

        height_h, height_f, p_h, p_f = calculate_heights_and_pressures(
            p_surf, height_surf, a_coeffs_137, b_coeffs_137, t_field, q_field,
        )
        ds_from_era5["height_h"][time_index] = height_h
        ds_from_era5["height_f"][time_index] = height_f
        ds_from_era5["p_h"][time_index] = p_h
        ds_from_era5["p_f"][time_index] = p_f


def era5_on_height_levels(ds_model_levels, heights_array):
    """Converts ERA5 model level data to data on height levels
    using Steffen interpolation"""
    if isinstance(heights_array[0], numbers.Integral):
        raise Exception("Heights need to be floating numbers, rather than integers")
    heights_coord = {
        "lev": ("lev", heights_array, {"long_name": "altitude", "units": "metres"},)
    }
    ds_height_levels = xr.Dataset(
        coords={
            "time": ds_model_levels.time,
            **heights_coord,
            "latitude": ds_model_levels.latitude,
            "longitude": ds_model_levels.longitude,
        }
    )
    time_steps = len(ds_model_levels["height_f"])
    shape_p_levels = np.shape(ds_model_levels["height_f"])
    shape_h_levels = (shape_p_levels[0],) + (len(heights_array),) + shape_p_levels[2:]
    for variable in ds_model_levels.variables:
        if ds_model_levels[variable].dims == (
            "time",
            "level",
            "latitude",
            "longitude",
        ):
            ds_height_levels[variable] = (
                ("time", "lev", "latitude", "longitude"),
                np.empty(shape_h_levels),
                ds_model_levels[variable].attrs,
            )
        elif "level" not in ds_model_levels[variable].dims:
            ds_height_levels[variable] = (
                ds_model_levels[variable].dims,
                ds_model_levels[variable],
                ds_model_levels[variable].attrs,
            )
    for time_index in range(time_steps):
        h_f_inverse = ds_model_levels["height_f"][time_index, ::-1, :, :].values
        h_h_inverse = ds_model_levels["height_h"][time_index, ::-1, :, :].values
        sea_mask = (
            (ds_model_levels["height_h"][time_index, -1, :, :].values < 5.0)
            * (ds_model_levels["height_h"][time_index, -1, :, :].values > 1.0e-6)
            * (ds_model_levels["lsm"][time_index, :, :].values < 0.2)
        )
        lower_extrapolation_with_mask = xr.where(
            sea_mask, -1.0e-6, ds_model_levels["height_h"][time_index, -1, :, :]
        ).values
        for variable in ds_model_levels.variables:
            if np.shape(ds_model_levels[variable]) == shape_p_levels:
                if variable in ["height_h", "p_h"]:
                    h_inverse = h_h_inverse
                else:
                    h_inverse = h_f_inverse
                field_p_levels = ds_model_levels[variable][
                    time_index, ::-1, :, :
                ].values
                if variable in ["p_h", "p_f", "height_h", "height_f"]:
                    ds_height_levels[variable][time_index] = steffen_3d(
                        field_p_levels,
                        h_inverse,
                        heights_array,
                        lower_extrapolation_with_mask,
                        lower_extrapolation_with_gradient=True,
                    )
                else:
                    ds_height_levels[variable][time_index] = steffen_3d(
                        field_p_levels,
                        h_inverse,
                        heights_array,
                        lower_extrapolation_with_mask,
                    )
    return ds_height_levels
