"""
ERA5 utilities that
- Add heights and pressures to an input data array on model levels
- Interpolate from model levels to constant heights
- Calculate gradients
- Subset domain

TODO
- Check code for surface height different from zero over ocean
- Check use of float vs double
- Move some functionality to more generic utilities
- Extend float 
- Find out how to get more efficient code for lower level variables
"""

import os
import numpy as np
import pandas as pd
import xarray as xr

# Optional numba dependency
try:
    from numba import njit

    print("Running with numba")
except ImportError:

    def njit(numba_function):
        """Dummy numba function"""
        return numba_function

    print("Running without numba")

# ECMWF CONSTANTS
rd = 287.06
rg = 9.80665
rv_over_rd_minus_one = 0.609133

# OTHER CONSTANTS NEED TO GO ELSEWHERE EVENTUALLY?
p_ref = 1.0e5
cp = 1004.0
rd_over_cp = rd / cp
p_ref_inv = 1.0 / p_ref
r_earth = 6371000.0

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


@njit
def steffen_3d(input_data, input_levels, output_level_array):
    """
    Performs Steffen interpolation on each individual column.
    Steffen, M. (1990). A simple method for monotonic interpolation
    in one dimension. Astronomy and Astrophysics, 239, 443.
    """
    k_max = input_data.shape[0]
    j_max = input_data.shape[1]
    i_max = input_data.shape[2]
    k_max_output = output_level_array.shape[0]
    k_max_minus = k_max - 1
    linear_slope = np.empty((k_max, j_max, i_max))
    output_data = np.empty((k_max_output, j_max, i_max))
    for i in range(i_max):
        for j in range(j_max):
            # first point
            delta_lower = input_levels[1, j, i] - input_levels[0, j, i]
            delta_upper = input_levels[2, j, i] - input_levels[1, j, i]
            if delta_lower < 0:
                raise Exception("Non-montonic increase in input_levels")
            if delta_upper < 0:
                raise Exception("Non-montonic increase in input_levels")
            slope_lower = (input_data[1, j, i] - input_data[0, j, i]) / delta_lower
            slope_upper = (input_data[2, j, i] - input_data[1, j, i]) / delta_upper

            weighted_slope = slope_lower * (
                1 + delta_lower / (delta_lower + delta_upper)
            ) - slope_upper * delta_lower / (delta_lower + delta_upper)
            if weighted_slope * slope_lower <= 0.0:
                linear_slope[0, j, i] = 0.0
            elif np.abs(weighted_slope) > 2 * np.abs(slope_lower):
                linear_slope[0, j, i] = 2.0 * slope_lower
            else:
                linear_slope[0, j, i] = weighted_slope

            # intermediate points
            for k in range(1, k_max_minus):
                delta_lower = input_levels[k, j, i] - input_levels[k - 1, j, i]
                delta_upper = input_levels[k + 1, j, i] - input_levels[k, j, i]
                slope_lower = (
                    input_data[k, j, i] - input_data[k - 1, j, i]
                ) / delta_lower
                slope_upper = (
                    input_data[k + 1, j, i] - input_data[k, j, i]
                ) / delta_upper
                weighted_slope = (
                    slope_lower * delta_upper + slope_upper * delta_lower
                ) / (delta_lower + delta_upper)

                if slope_lower * slope_upper <= 0.0:
                    linear_slope[k, j, i] = 0.0
                elif np.abs(weighted_slope) > 2.0 * np.abs(slope_lower):
                    linear_slope[k, j, i] = np.copysign(2.0, slope_lower) * min(
                        np.abs(slope_lower), np.abs(slope_upper)
                    )
                elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
                    linear_slope[k, j, i] = np.copysign(2.0, slope_lower) * min(
                        np.abs(slope_lower), np.abs(slope_upper)
                    )
                else:
                    linear_slope[k, j, i] = weighted_slope

            # last point
            delta_lower = (
                input_levels[k_max_minus - 1, j, i]
                - input_levels[k_max_minus - 2, j, i]
            )
            delta_upper = (
                input_levels[k_max_minus, j, i] - input_levels[k_max_minus - 1, j, i]
            )
            slope_lower = (
                input_data[k_max_minus - 1, j, i] - input_data[k_max_minus - 2, j, i]
            ) / delta_lower
            slope_upper = (
                input_data[k_max_minus, j, i] - input_data[k_max_minus - 1, j, i]
            ) / delta_upper
            weighted_slope = slope_upper * (
                1 + delta_upper / (delta_upper + delta_lower)
            ) - slope_lower * delta_upper / (delta_upper + delta_lower)
            if weighted_slope * slope_upper <= 0.0:
                linear_slope[k_max_minus, j, i] = 0.0
            elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
                linear_slope[k_max_minus, j, i] = 2.0 * slope_upper
            else:
                linear_slope[k_max_minus, j, i] = weighted_slope

            # loop over output points
            k_temp = 0
            for k_out in range(k_max_output):
                while (k_temp < k_max) and (
                    input_levels[k_temp, j, i] < output_level_array[k_out]
                ):
                    k_temp = k_temp + 1
                if k_temp > 0 and k_temp < k_max:
                    k_high = k_temp
                    k_low = k_high - 1
                    delta = input_levels[k_high, j, i] - input_levels[k_low, j, i]
                    slope = (input_data[k_high, j, i] - input_data[k_low, j, i]) / delta
                    a = (
                        linear_slope[k_low, j, i]
                        + linear_slope[k_high, j, i]
                        - 2 * slope
                    ) / (delta * delta)
                    b = (
                        3 * slope
                        - 2 * linear_slope[k_low, j, i]
                        - linear_slope[k_high, j, i]
                    ) / delta
                    c = linear_slope[k_low, j, i]
                    d = input_data[k_low, j, i]
                    t = output_level_array[k_out] - input_levels[k_low, j, i]
                    t_2 = t * t
                    t_3 = t_2 * t
                    output_data[k_out, j, i] = a * t_3 + b * t_2 + c * t + d
                else:
                    output_data[k_out, j, i] = np.nan
    return output_data


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
    return ds_from_era5


def era5_on_height_levels(ds_pressure_levels, heights_array):
    """Converts ERA5 model level data to data on height levels
    using Steffen interpolation"""
    heights_coord = {
        "height": (
            "height",
            heights_array,
            {"long_name": "height above sea level", "units": "metres"},
        )
    }
    ds_height_levels = xr.Dataset(
        coords={
            "time": ds_pressure_levels.time,
            **heights_coord,
            "latitude": ds_pressure_levels.latitude,
            "longitude": ds_pressure_levels.longitude,
        }
    )
    time_steps = len(ds_pressure_levels["height_f"])
    shape_p_levels = np.shape(ds_pressure_levels["height_f"])
    shape_h_levels = (shape_p_levels[0],) + (len(heights_array),) + shape_p_levels[2:]
    shape_no_levels = (shape_p_levels[0],) + shape_p_levels[2:]
    for variable in ds_pressure_levels.variables:
        if ds_pressure_levels[variable].dims == (
            "time",
            "level",
            "latitude",
            "longitude",
        ):
            ds_height_levels[variable] = (
                ("time", "height", "latitude", "longitude"),
                np.empty(shape_h_levels),
                ds_pressure_levels[variable].attrs,
            )
        elif "level" not in ds_pressure_levels[variable].dims:
            ds_height_levels[variable] = (
                ds_pressure_levels[variable].dims,
                ds_pressure_levels[variable],
                ds_pressure_levels[variable].attrs,
            )
    for time_index in range(time_steps):
        h_f_inverse = ds_pressure_levels["height_f"][time_index, ::-1, :, :].values
        h_h_inverse = ds_pressure_levels["height_h"][time_index, ::-1, :, :].values
        for variable in ds_pressure_levels.variables:
            if np.shape(ds_pressure_levels[variable]) == shape_p_levels:
                if variable in ["height_h", "p_h"]:
                    h_inverse = h_h_inverse
                else:
                    h_inverse = h_f_inverse
                field_p_levels = ds_pressure_levels[variable][
                    time_index, ::-1, :, :
                ].values
                ds_height_levels[variable][time_index] = steffen_3d(
                    field_p_levels, h_inverse, heights_array
                )
    return ds_height_levels

def add_lower_level_variables(ds_pressure_levels):
    """Extend below the lowest full level to the corresponding
    half-level. Over sea, mak sure to extend to below zero.
    Unfortunately, 'normal concatenation' of data-arrays struggles
    perform well. The current code seems pretty memory-inefficient"""
    sea_mask = (
        (ds_pressure_levels["height_h"][:, -1, :, :].values < 5.0)
        * (ds_pressure_levels["height_h"][:, -1, :, :].values > 1.0e-6)
        * (ds_pressure_levels["lsm"].values < 0.2)
    )
    level_data_dims = (
        "time",
        "level",
        "latitude",
        "longitude",
    )
    levels_extended = {
        "level": (
            "level",
            np.concatenate(
                (
                    ds_pressure_levels["level"].values,
                    [ds_pressure_levels["level"].values[-1] + 1],
                )
            ).astype("int32"),
            ds_pressure_levels["level"].attrs,
        )
    }
    coords_extended = {
        "time": ds_pressure_levels.time,
        **levels_extended,
        "latitude": ds_pressure_levels.latitude,
        "longitude": ds_pressure_levels.longitude,
    }
    ds_extended = xr.Dataset(coords=coords_extended)
    for variable in ds_pressure_levels.variables:
        if ds_pressure_levels[variable].dims == level_data_dims:
            ds_extended[variable] = (
                coords_extended,
                np.concatenate(
                    (
                        ds_pressure_levels[variable][:, :, :, :],
                        ds_pressure_levels[variable][:, [-1], :, :],
                    ),
                    axis=1,
                ),
            )
        elif "level" not in ds_pressure_levels[variable].dims:
            ds_extended[variable] = (
                ds_pressure_levels[variable].dims,
                ds_pressure_levels[variable],
                ds_pressure_levels[variable].attrs,
            )
    ds_extended["height_f"][:, -1, :, :] = ds_extended["height_h"][:, -1, :, :]
    ds_extended["p_f"][:, -1, :, :] = ds_extended["p_h"][:, -1, :, :]
    ds_extended["height_h"][:, -1, :, :] = ds_extended["height_h"][:, -2, :, :] - (
        ds_extended["height_h"][:, -3, :, :] - ds_extended["height_h"][:, -2, :, :]
    )
    ds_extended["p_h"][:, -1, :, :] = ds_extended["p_h"][:, -2, :, :] - (
        ds_extended["p_h"][:, -3, :, :] - ds_extended["p_h"][:, -2, :, :]
    )
    ds_extended["height_f"][:, -1, :, :] = xr.where(
        sea_mask, -1.0e-6, ds_extended["height_f"][:, -1, :, :]
    )
    slope = (ds_extended["p_h"][:, -3, :, :] - ds_extended["p_h"][:, -2, :, :]) / (
        ds_extended["height_h"][:, -3, :, :] - ds_extended["height_h"][:, -2, :, :]
    )
    ds_extended["p_f"][:, -1, :, :] = xr.where(
        sea_mask,
        ds_extended["p_h"][:, -2, :, :]
        - (ds_extended["height_h"][:, -2, :, :] - ds_extended["height_f"][:, -1, :, :])
        * slope,
        ds_extended["p_f"][:, -1, :, :],
    )
    return ds_extended


def add_auxiliary_variable(ds_level_2, var):
    """Alternative: equations could be separated out to utility"""
    if var == "theta":
        attr_dict = {"units": "K", "long_name": "potential temperature"}
        ds_level_2[var] = (
            ds_level_2["t"] * (ds_level_2["p_f"] * p_ref_inv) ** rd_over_cp
        )
    else:
        raise NotImplementedError("Variable not implemented")
    ds_level_2[var] = ds_level_2[var].assign_attrs(**attr_dict)
    return ds_level_2


def add_auxiliary_variables(ds_level_1, list_of_vars):
    """Add variables defined in list to dictionary"""
    for var in list_of_vars:
        add_auxiliary_variable(ds_level_1, var)
    return ds_level_1


def era_5_normalise_longitude(ds_to_normalise):
    """Normalise lognitudes"""
    ds_to_normalise.coords["longitude"] = np.round(
        ds_to_normalise.coords["longitude"] % 360.0, decimals=4
    )
    return ds_to_normalise


def era_5_subset(ds_full, dictionary):
    """Utility to select era5 data by latitude and longitude
    Note: data order is North to South"""
    ds_subset = ds_full.sel(
        latitude=slice(dictionary["lat_max"], dictionary["lat_min"]),
        longitude=slice(dictionary["lon_min"] % 360, dictionary["lon_max"] % 360),
    )
    return ds_subset

def era5_single_point(ds_domain, dictionary):
    """for a single variable, extract a local profile"""
    ds_at_location = ds_domain.sel(
        latitude=dictionary["lat"], longitude=dictionary["lon"] % 360, method="nearest"
    )
    return ds_at_location

def era5_mask(ds_to_mask, dictionary):
    # Only use ocean points, ensure it can be used after before or after array extensions
    if dictionary["mask"] == "ocean":
        mask = (ds_to_mask["z"] < 5.0 * rg) * (ds_to_mask["lsm"].values < 0.2)
    return mask

def era5_weighted(ds_to_weigh, dictionary):
    """adds weights to dictionary"""
    if "weights" in dictionary:
        if dictionary["weights"] == "area":
            ds_to_weigh.weigths = np.cos(np.deg2rad(ds_to_weigh.lat_meshgrid))
        else:
            raise Exception("weight strategy not implemented")
    return ds_to_weigh


def era5_box_mean(ds_box, dictionary):
    """
    - Only use columns where the first level is higher than the local first level
    in the location of interest
    - Option to weight by box size?
    """
    era5_weighted(ds_box, dictionary)
    if "mask" in dictionary:
        mask = era5_mask(ds_box, dictionary)
        ds_mean = ds_box.where(mask).mean(("latitude", "longitude"))
    else:
        ds_mean = ds_box.mean(("latitude", "longitude"))
    return ds_mean

def era5_add_lat_lon_meshgrid(ds_to_extend):
    """add a lat, lon meshgrid to a dataset, useful for gradients with filter"""
    lon_mesh, lat_mesh = np.meshgrid(ds_to_extend.longitude, ds_to_extend.latitude)
    ds_to_extend["lon_meshgrid"] = (
        ("latitude", "longitude"),
        lon_mesh,
        ds_to_extend["longitude"].attrs,
    )
    ds_to_extend["lat_meshgrid"] = (
        ("latitude", "longitude"),
        lat_mesh,
        ds_to_extend["latitude"].attrs,
    )
    return ds_to_extend

def dist(ds_1, ds_2):
    pi = np.pi
    lon1 = ds_1["lon_meshgrid"].values * (2 * pi / 360)
    lon2 = ds_2["lon_meshgrid"].values * (2 * pi / 360)
    dlon = lon2 - lon1
    lat1 = ds_1["lat_meshgrid"].values * (2 * pi / 360)
    lat2 = ds_2["lat_meshgrid"].values * (2 * pi / 360)
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    dist = r_earth * c
    return dist

def era5_boundary_gradients(ds_box, variable, dictionary):
    """ Calculate gradients using haversine function
    How to deal with missing data?
    Weight by box size?
    Gradients from boundary values"""
    # use argmin instead ??
    # ~ left = left.where(np.nanmin(ds_box["height_f"][:, :, :, :]) < 20.0)
    # ~ left = left.where(left.longitude > [dictionary["lon_min"] % 360])
    # left = left.sel(latitude=slice(dictionary["lat_max"], dictionary["lat_min"]))
    left = ds_box.min("longitude", skipna=True)
    right = ds_box.max("longitude", skipna=True)
    top = ds_box.max("latitude", skipna=True)
    bottom = ds_box.min("latitude", skipna=True)
    x_gradient = np.mean(
        (right[variable].values - left[variable].values) / dist(left, right), axis=2
    )
    y_gradient = np.mean(
        (top[variable].values - bottom[variable].values) / dist(top, bottom), axis=2
    )
    return x_gradient, y_gradient


def era5_regression_gradients(ds_box, variable, dictionary):
    """ Calculate gradients using haversine function
    From regression, using the normal equation"""
    pi = np.pi
    lon1 = dictionary["lon"] * (2 * pi / 360)
    lon2 = ds_box["lon_meshgrid"].values * (2 * pi / 360)
    dlon = lon2 - lon1
    lat1 = dictionary["lat"] * (2 * pi / 360)
    lat2 = ds_box["lat_meshgrid"].values * (2 * pi / 360)
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    dist = r_earth * c
    theta = np.arctan2(
        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon),
        np.sin(dlon) * np.cos(lat2),
    )
    x = dist * np.cos(theta)
    y = dist * np.sin(theta)
    x_flat = x.flatten()
    y_flat = y.flatten()
    ones_flat = np.ones(np.shape(x_flat))
    len_temp = np.shape(ds_box[variable])[0]
    len_levels = np.shape(ds_box[variable])[1]
    x_gradient_array = np.empty((len_temp, len_levels))
    y_gradient_array = np.empty((len_temp, len_levels))
    for this_time in range(len_temp):
        for this_level in range(len_levels):
            data_flat = ds_box[variable][this_time, this_level, :, :].values.flatten()
            data_flat_filter = data_flat[~np.isnan(data_flat)][:,None]
            x_flat_filter = x_flat[~np.isnan(data_flat)][:,None]
            y_flat_filter = y_flat[~np.isnan(data_flat)][:,None]
            ones_flat_filter= ones_flat[~np.isnan(data_flat)][:,None]
            oxy_mat = np.hstack((ones_flat_filter, x_flat_filter, y_flat_filter))
            theta = np.dot(
                np.dot(np.linalg.pinv(np.dot(oxy_mat.transpose(), oxy_mat)), oxy_mat.transpose()),
                data_flat_filter,
            )
            x_gradient_array[this_time, this_level] = theta[1]
            y_gradient_array[this_time, this_level] = theta[2]
    return x_gradient_array, y_gradient_array

def era5_add_gradients_to_variables(ds_level_1, list_of_vars, dictionary):
    """Add variables defined in list to dictionary"""
    ds_out = xr.Dataset(coords={"time": ds_level_1.time, "height": ds_level_1.height})
    for variable in list_of_vars:
        if dictionary["gradients_strategy"] in ["regression", "both"]:
            x_gradient_array, y_gradient_array = era5_regression_gradients(
                ds_level_1, variable, dictionary
            )
        elif dictionary["gradients_strategy"] == "boundary":
            x_gradient_array, y_gradient_array = era5_boundary_gradients(
                ds_level_1, variable, dictionary
            )
        else:
            raise NotImplementedError("Gradients strategy not implemented")
        ds_out["d" + variable + "dx"] = (
            ("time", "height"),
            x_gradient_array,
        )
        ds_out["d" + variable + "dy"] = (
            ("time", "height"),
            y_gradient_array,
        )
        if dictionary["gradients_strategy"] == "both":
            x_gradient_array, y_gradient_array = era5_boundary_gradients(
                ds_level_1, variable, dictionary
            )
            ds_out["d" + variable + "dx_bound"] = (
                ("time", "height"),
                x_gradient_array,
            )
            ds_out["d" + variable + "dy_bound"] = (
                ("time", "height"),
                y_gradient_array,
            )
    return ds_out

if __name__ == "__main__":
    # just building an example here, with lazy loading
    files_model_an = "output_domains/model_an_*_eurec4a_circle_eul_domain.nc"
    files_single_an = "output_domains/single_an_*_eurec4a_circle_eul_domain.nc"
    files_model_fc = "output_domains/model_fc_*_eurec4a_circle_eul_domain.nc"
    files_single_fc = "output_domains/single_fc_*_eurec4a_circle_eul_domain.nc"
    ds_model_an = xr.open_mfdataset(files_model_an, combine="by_coords")
    ds_model_an = ds_model_an.drop_vars("z")
    ds_single_an = xr.open_mfdataset(files_single_an, combine="by_coords")
    ds_model_fc = xr.open_mfdataset(files_model_fc, combine="by_coords")
    ds_single_fc = xr.open_mfdataset(files_single_fc, combine="by_coords")
    ds_list = [ds_model_an, ds_single_an, ds_model_fc, ds_single_fc]
    for this_ds in ds_list:
        era_5_normalise_longitude(this_ds)
    ds_merged = xr.merge(ds_list)
    ds_time = ds_merged.isel(time=[24])
    lats_lons_dictionary = {
        "lat_min": 11.3,
        "lat_max": 15.3,
        "lon_min": -59.717,
        "lon_max": -55.717,
        "lat": 13.3,
        "lon": -57.717,
        "gradients_strategy": "both",
    }
    ds_smaller = era_5_subset(ds_time, lats_lons_dictionary)
    add_heights_and_pressures(ds_smaller)
    ds_extended = add_lower_level_variables(ds_smaller)
    add_auxiliary_variables(ds_extended, ["theta"])
    heights_array = np.arange(0, 5000, 40)
    ds_time_height = era5_on_height_levels(ds_extended, heights_array)
    era5_add_lat_lon_meshgrid(ds_time_height)
    ds_gradients = era5_add_gradients_to_variables(
        ds_time_height, ["u", "v", "p_f", "theta"], lats_lons_dictionary
    )
    ds_gradients.to_netcdf("ds_gradients.nc")
    # era5_regression_gradient(ds_time_height, "u", lats_lons_dictionary)
    # era5_box_average_gradients(ds_time, lats_lons_dictionary)
    # ~ ds_era5_single_point = era5_single_point(ds_time_height, lats_lons_dictionary)
    # ~ ds_era5_single_point.to_netcdf("ds_profile.nc")
    # ~ ds_era5_mean = era5_box_mean(ds_time_height, lats_lons_dictionary)
    # ~ ds_era5_mean.to_netcdf("ds_mean.nc")
