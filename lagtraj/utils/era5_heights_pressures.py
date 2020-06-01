"""
ERA5 utilities that can
- Add heights and pressures to an input data array on model levels
- Interpolate from model levels to constant height levels (using Steffen interpolation)
- Calculate gradients using boundary values or regression method
- Extract local profiles and mean profiles
- Subselect a domain
- Filter/mask data: e.g. "ocean values only"
- Add auxiliary variables

TODO
- Optimise code (note that interpolation is currently expensive, possibly because coordinates are not assumed to be ordered)
- Fix potential nan-gradients?
- Add more auxiliary variables
- Move some functionality (e.g. auxiliary variables) to more generic utilities
- Test/develop way of dealing with -180 degrees.
  Note: HIGH-TUNE/DEPHY prefers longitude between -180..180?
- Implement HIGH-TUNE conventions, variable renaming, attributes
- HIGH-TUNE/DEPHY needs netcdf3?
- Discuss need to check/convert float vs double (HIGH-TUNE/DEPHY expects double)
- Discuss data filter and weight procedures
- Use more exact means and gradients based on interpolation/weights?
- Compare regression and boundary gradients for data
- Look into other mean/gradient techniques (e.g. Gaussian weighted, see CSET code).
- Further documentation
- Keep checking against cf conventions
"""

import os
import numbers
import datetime
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
pi = np.pi
Omega = 7.2921150e-5

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
def steffen_3d(
    input_data,
    input_levels,
    output_level_array,
    lower_extrapolation_surface,
    lower_extrapolation_with_gradient=False,
):
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
                if 0 < k_temp < k_max:
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
                    t_1 = output_level_array[k_out] - input_levels[k_low, j, i]
                    t_2 = t_1 * t_1
                    t_3 = t_2 * t_1
                    output_data[k_out, j, i] = a * t_3 + b * t_2 + c * t_1 + d
                elif (k_temp == 0) and (
                    output_level_array[k_out] >= lower_extrapolation_surface[j, i]
                ):
                    if lower_extrapolation_with_gradient:
                        output_data[k_out, j, i] = input_data[0, j, i] + linear_slope[
                            0, j, i
                        ] * (output_level_array[k_out] - input_levels[0, j, i])
                    else:
                        output_data[k_out, j, i] = input_data[0, j, i]
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


def era5_on_height_levels(ds_pressure_levels, heights_array):
    """Converts ERA5 model level data to data on height levels
    using Steffen interpolation"""
    if isinstance(heights_array[0], numbers.Integral):
        raise Exception("Heights need to be floating numbers, rather than integers")
    heights_coord = {
        "lev": ("lev", heights_array, {"long_name": "altitude", "units": "metres"},)
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
    for variable in ds_pressure_levels.variables:
        if ds_pressure_levels[variable].dims == (
            "time",
            "level",
            "latitude",
            "longitude",
        ):
            ds_height_levels[variable] = (
                ("time", "lev", "latitude", "longitude"),
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
        sea_mask = (
            (ds_pressure_levels["height_h"][time_index, -1, :, :].values < 5.0)
            * (ds_pressure_levels["height_h"][time_index, -1, :, :].values > 1.0e-6)
            * (ds_pressure_levels["lsm"][time_index, :, :].values < 0.2)
        )
        lower_extrapolation_with_mask = xr.where(
            sea_mask, -1.0e-6, ds_pressure_levels["height_h"][time_index, -1, :, :]
        ).values
        for variable in ds_pressure_levels.variables:
            if np.shape(ds_pressure_levels[variable]) == shape_p_levels:
                if variable in ["height_h", "p_h"]:
                    h_inverse = h_h_inverse
                else:
                    h_inverse = h_f_inverse
                field_p_levels = ds_pressure_levels[variable][
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


def add_auxiliary_variable(ds_to_expand, var):
    """Adds auxiliary variables to arrays.
    Alternatively, the equations could be separated out to another utility
    I think this may be adding a 'black box layer' though
    To be discussed"""
    if var == "theta":
        attr_dict = {"units": "K", "long_name": "potential temperature"}
        ds_to_expand[var] = (
            ds_to_expand["t"] * (ds_to_expand["p_f"] * p_ref_inv) ** rd_over_cp
        )
    elif var == "rho":
        attr_dict = {"units": "kg m**-3", "long_name": "density"}
        ds_to_expand[var] = ds_to_expand["p_f"] / (
            rd * ds_to_expand["t"] * (1.0 + rv_over_rd_minus_one * ds_to_expand["q"])
        )
    else:
        raise NotImplementedError("Variable not implemented")
    ds_to_expand[var] = ds_to_expand[var].assign_attrs(**attr_dict)


def add_auxiliary_variables(ds_to_expand, list_of_vars):
    """Wrapper for auxiliary variable calculation"""
    for var in list_of_vars:
        add_auxiliary_variable(ds_to_expand, var)


def longitude_set_meridian(longitude):
    """Sets longitude to be between -180 and 180 degrees"""
    return (longitude + 180.0) % 360.0 - 180.0


def era_5_normalise_longitude(ds_to_normalise):
    """Normalise longitudes to be between 0 and 360 degrees
    This is needed because these are stored differently in the surface
    and model level data. Rounding up to 4 decimals seems to work for now,
    with more decimals misalignment has happenend. Would be good to sort
    out why this is the case.
    """
    ds_to_normalise.coords["longitude"] = (
        "longitude",
        np.round(
            longitude_set_meridian(ds_to_normalise.coords["longitude"]), decimals=4
        ),
        ds_to_normalise.coords["longitude"].attrs,
    )
    return ds_to_normalise


def era_5_subset(ds_full, dictionary):
    """Utility to select era5 data by latitude and longitude
    Note: data order is North to South"""
    ds_subset = ds_full.sel(
        latitude=slice(dictionary["lat_max"], dictionary["lat_min"]),
        longitude=slice(
            longitude_set_meridian(dictionary["lon_min"]),
            longitude_set_meridian(dictionary["lon_max"]),
        ),
    )
    return ds_subset


def era5_single_point(ds_domain, dictionary):
    """Extracts a local profile at the nearest point"""
    ds_at_location = ds_domain.sel(
        latitude=dictionary["lat"],
        longitude=longitude_set_meridian(dictionary["lon"]),
        method="nearest",
    )
    return ds_at_location


def era5_interp_column(ds_domain, lat_to_interp, lon_to_interp):
    """Returns the dataset interpolated to given latitude and longitude
    with latitude and longitude dimensions retained"""
    ds_at_location = ds_domain.interp(
        latitude=[lat_to_interp], longitude=[longitude_set_meridian(lon_to_interp)]
    )
    return ds_at_location


def era5_mask(ds_to_mask, dictionary):
    """Returns a lat-lon mask"""
    # Only use ocean points, ensure it can be used after before or after array extensions
    if dictionary["mask"] == "ocean":
        mask = (ds_to_mask["z"] < 5.0 * rg) * (ds_to_mask["lsm"].values < 0.2)
    return mask


def era5_weighted(ds_to_weigh, dictionary):
    """Adds weights to dictionary"""
    if "weights" in dictionary:
        if dictionary["weights"] == "area":
            ds_to_weigh.weigths = np.cos(np.deg2rad(ds_to_weigh.lat_meshgrid))
        else:
            raise Exception("weight strategy not implemented")


def era5_box_mean(ds_box, dictionary):
    """
    Calculates mean over a data_set.
    - Only use columns where the first level is higher than the local first level
    in the location of interest
    - Option to weight by box size?
    """
    era5_weighted(ds_box, dictionary)
    if "mask" in dictionary:
        mask = era5_mask(ds_box, dictionary)
        ds_mean = ds_box.where(mask).mean(("latitude", "longitude"), keep_attrs=True)
    else:
        ds_mean = ds_box.mean(("latitude", "longitude"), keep_attrs=True)
    return ds_mean


def era5_add_lat_lon_meshgrid(ds_to_extend):
    """Adds a [lat, lon] meshgrid to a dataset, useful for gradients in order to work
    around nan values on edge"""
    lon_mesh, lat_mesh = np.meshgrid(ds_to_extend.longitude, ds_to_extend.latitude)
    ds_to_extend["lon_meshgrid"] = (
        ("latitude", "longitude"),
        lon_mesh,
        {
            "long_name": "longitude on meshgrid",
            "units": ds_to_extend["longitude"].units,
        },
    )
    ds_to_extend["lat_meshgrid"] = (
        ("latitude", "longitude"),
        lat_mesh,
        {"long_name": "latitude on meshgrid", "units": ds_to_extend["latitude"].units},
    )


def calc_haver_dist(lat1, lon1, lat2, lon2):
    """Calculates distance given pairs of latitude and longitude in radians"""
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    haver = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    arc_dist = 2 * np.arctan2(np.sqrt(haver), np.sqrt(1.0 - haver))
    haver_dist = r_earth * arc_dist
    return haver_dist


def calc_lat_lon_angle(lat1, lon1, lat2, lon2):
    """Calculates angle given pairs of latitude and longitude in radians"""
    dlon = lon2 - lon1
    lat_lon_angle = np.arctan2(
        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon),
        np.sin(dlon) * np.cos(lat2),
    )
    return lat_lon_angle


def lat_dist(lat1, lat2):
    """Distance for difference in latitude only, at constant longitude"""
    return np.deg2rad(lat2 - lat1) * r_earth


def lon_dist(lon1, lon2, lat):
    """Distance for difference in longitude, for given latitude"""
    dlon_rad = np.deg2rad(lon2 - lon1)
    dlon_rad = (dlon_rad + pi) % (2 * pi) - pi
    return np.cos(np.deg2rad(lat)) * dlon_rad * r_earth


def lon_dist_from_meshgrids(ds_1, ds_2):
    """Calculates distances between two datasets of the same shape"""
    lat1_mg = ds_1["lat_meshgrid"].values
    lon1_mg = ds_1["lon_meshgrid"].values
    lat2_mg = ds_2["lat_meshgrid"].values
    lon2_mg = ds_2["lon_meshgrid"].values
    return lon_dist(lon1_mg, lon2_mg, 0.5 * (lat1_mg + lat2_mg)).flatten()


def lat_dist_from_meshgrids(ds_1, ds_2):
    """Calculates distances between two datasets of the same shape"""
    lat1_mg = ds_1["lat_meshgrid"].values
    lat2_mg = ds_2["lat_meshgrid"].values
    return lat_dist(lat1_mg, lat2_mg).flatten()


@njit
def boundary_gradients(x_array, y_array, val_array):
    """Numba function to calculate gradients, dismissing filtered points"""
    len_temp = np.shape(val_array)[0]
    len_levels = np.shape(val_array)[1]
    len_lats = np.shape(val_array)[2]
    len_lons = np.shape(val_array)[3]
    x_gradient_array = np.empty((len_temp, len_levels))
    y_gradient_array = np.empty((len_temp, len_levels))
    dval_dx = np.empty((len_lats))
    dval_dy = np.empty((len_lons))
    for this_time in range(len_temp):
        for this_level in range(len_levels):
            for this_lat in range(len_lats):
                vals = val_array[this_time, this_level, this_lat, :].flatten()
                x_at_lat = x_array[this_lat, :].flatten()
                vals_filtered = vals[~np.isnan(vals)]
                x_filtered = x_at_lat[~np.isnan(vals)]
                dvals = vals_filtered[-1] - vals_filtered[0]
                dval_dx[this_lat] = dvals / (x_filtered[-1] - x_filtered[0])
            for this_lon in range(len_lons):
                vals = val_array[this_time, this_level, :, this_lon].flatten()
                y_at_lat = y_array[:, this_lon].flatten()
                vals_filtered = vals[~np.isnan(vals)]
                y_filtered = y_at_lat[~np.isnan(vals)]
                dvals = vals_filtered[-1] - vals_filtered[0]
                dval_dy[this_lon] = dvals / (y_filtered[-1] - y_filtered[0])
            x_gradient_array[this_time, this_level] = np.mean(dval_dx)
            y_gradient_array[this_time, this_level] = np.mean(dval_dy)
    return x_gradient_array, y_gradient_array


def era5_boundary_gradients(ds_box, variable, dictionary):
    """ Calculate gradients from boundary values
    Using distances along latitude and longitude axes
    Weight by box size?"""
    ds_filtered = ds_box
    if "mask" in dictionary:
        mask = era5_mask(ds_box, dictionary)
        ds_filtered = ds_filtered.where(mask)
    lat1_point = dictionary["lat"]
    lon1_point = longitude_set_meridian(dictionary["lon"])
    lat2_mg = ds_filtered["lat_meshgrid"].values
    lon2_mg = ds_filtered["lon_meshgrid"].values
    x_array = lon_dist(lon1_point, lon2_mg, lat2_mg)
    y_array = lat_dist(lat1_point, lat2_mg)
    val_array = ds_filtered[variable].values
    return boundary_gradients(x_array, y_array, val_array)


def era5_regression_gradients(ds_box, variable, dictionary):
    """ Calculate gradients function using local coordinate system and
    Regression, using the normal equation"""
    ds_filtered = ds_box
    if "mask" in dictionary:
        mask = era5_mask(ds_box, dictionary)
        ds_filtered = ds_filtered.where(mask)
    lat1_point = dictionary["lat"]
    lon1_point = longitude_set_meridian(dictionary["lon"])
    lat2_mg = ds_filtered["lat_meshgrid"].values
    lon2_mg = ds_filtered["lon_meshgrid"].values
    x_array = lon_dist(lon1_point, lon2_mg, lat1_point)
    y_array = lat_dist(lat1_point, lat2_mg)
    x_flat = x_array.flatten()
    y_flat = y_array.flatten()
    ones_flat = np.ones(np.shape(x_flat))
    len_temp = np.shape(ds_box[variable])[0]
    len_levels = np.shape(ds_box[variable])[1]
    x_gradient_array = np.empty((len_temp, len_levels))
    y_gradient_array = np.empty((len_temp, len_levels))
    for this_time in range(len_temp):
        for this_level in range(len_levels):
            data_flat = ds_filtered[variable][
                this_time, this_level, :, :
            ].values.flatten()
            data_flat_filter = data_flat[~np.isnan(data_flat)][:, None]
            x_flat_filter = x_flat[~np.isnan(data_flat)][:, None]
            y_flat_filter = y_flat[~np.isnan(data_flat)][:, None]
            ones_flat_filter = ones_flat[~np.isnan(data_flat)][:, None]
            oxy_mat = np.hstack((ones_flat_filter, x_flat_filter, y_flat_filter))
            theta = np.dot(
                np.dot(
                    np.linalg.pinv(np.dot(oxy_mat.transpose(), oxy_mat)),
                    oxy_mat.transpose(),
                ),
                data_flat_filter,
            )
            x_gradient_array[this_time, this_level] = theta[1]
            y_gradient_array[this_time, this_level] = theta[2]
    return x_gradient_array, y_gradient_array


def era5_gradients(ds_field, list_of_vars, dictionary):
    """Add variables defined in list to dictionary"""
    ds_out = xr.Dataset(coords={"time": ds_field.time, "lev": ds_field.lev})
    for variable in list_of_vars:
        if dictionary["gradients_strategy"] in ["regression", "both"]:
            x_gradient_array, y_gradient_array = era5_regression_gradients(
                ds_field, variable, dictionary
            )
        elif dictionary["gradients_strategy"] == "boundary":
            x_gradient_array, y_gradient_array = era5_boundary_gradients(
                ds_field, variable, dictionary
            )
        else:
            raise NotImplementedError("Gradients strategy not implemented")
        ds_out["d" + variable + "dx"] = (
            ("time", "lev"),
            x_gradient_array,
            {
                "long_name": ds_field[variable].long_name + " x-gradient",
                "units": ds_field[variable].units + " m**-1",
            },
        )
        ds_out["d" + variable + "dy"] = (
            ("time", "lev"),
            y_gradient_array,
            {
                "long_name": ds_field[variable].long_name + " y-gradient",
                "units": ds_field[variable].units + " m**-1",
            },
        )
        if dictionary["gradients_strategy"] == "both":
            x_gradient_array, y_gradient_array = era5_boundary_gradients(
                ds_field, variable, dictionary
            )
            ds_out["d" + variable + "dx_bound"] = (
                ("time", "lev"),
                x_gradient_array,
                {
                    "long_name": ds_field[variable].long_name
                    + " x-gradient (boundaries)",
                    "units": ds_field[variable].units + " m**-1",
                },
            )
            ds_out["d" + variable + "dy_bound"] = (
                ("time", "lev"),
                y_gradient_array,
                {
                    "long_name": ds_field[variable].long_name
                    + " y-gradient (boundaries)",
                    "units": ds_field[variable].units + " m**-1",
                },
            )
    return ds_out


def era5_adv_tendencies(ds_profile, list_of_vars, dictionary):
    """Add variables defined in list to dictionary"""
    ds_out = xr.Dataset(coords={"time": ds_profile.time, "lev": ds_profile.lev})
    for variable in list_of_vars:
        tendency_array = (
            (ds_profile["u"].values - dictionary["u_traj"])
            * ds_profile["d" + variable + "dx"].values
            + (ds_profile["v"].values - dictionary["v_traj"])
            * ds_profile["d" + variable + "dy"].values
        )
        ds_out[variable + "_advtend"] = (
            ("time", "lev"),
            tendency_array,
            {
                "long_name": ds_profile[variable].long_name + " tendency (advection)",
                "units": ds_profile[variable].units + " s**-1",
            },
        )
        if dictionary["gradients_strategy"] == "both":
            tendency_array = (
                (ds_profile["u"].values - dictionary["u_traj"])
                * ds_profile["d" + variable + "dx_bound"].values
                + (ds_profile["v"].values - dictionary["v_traj"])
                * ds_profile["d" + variable + "dy_bound"].values
            )
            ds_out[variable + "_advtend_bound"] = (
                ("time", "lev"),
                tendency_array,
                {
                    "long_name": ds_profile[variable].long_name
                    + " tendency (advection, boundaries)",
                    "units": ds_profile[variable].units + " s**-1",
                },
            )

    return ds_out


def add_geowind_around_centre(ds_profile, dictionary):
    """Calculates the geostrophic wind at the centre of the box, based
    on mean profiles of density and gradients of pressure"""
    lat_centre = dictionary["lat"]
    f_cor = 2.0 * Omega * np.sin(np.deg2rad(lat_centre))
    u_geo = -(1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdy"]
    v_geo = (1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdx"]
    ds_profile["ug"] = (
        ("time", "lev"),
        u_geo,
        {"long_name": "U component of geostrophic wind", "units": "m s**-1",},
    )
    ds_profile["vg"] = (
        ("time", "lev"),
        v_geo,
        {"long_name": "V component of geostrophic wind", "units": "m s**-1",},
    )
    if dictionary["gradients_strategy"] == "both":
        u_geo_bound = -(1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdy_bound"]
        v_geo_bound = (1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdx_bound"]
        ds_profile["ug_bound"] = (
            ("time", "lev"),
            u_geo_bound,
            {
                "long_name": "U component of geostrophic wind (boundaries)",
                "units": "m s**-1",
            },
        )
        ds_profile["vg_bound"] = (
            ("time", "lev"),
            v_geo_bound,
            {
                "long_name": "V component of geostrophic wind (boundaries)",
                "units": "m s**-1",
            },
        )


def trace_one_way(lat, lon, u_traj, v_traj, d_time, lforward=True):
    """calculates previous position given lat,lon,u_traj,v_traj, and d_time
    explicitly set d_time positive, to prevent accidental backward trajectory"""
    if d_time < 0.0:
        raise Exception("Expecting positive d_time in back-tracing")
    theta = np.arctan2(v_traj, u_traj) % (2 * pi)
    if lforward:
        bearing = ((pi / 2) - theta) % (2 * pi)
    else:
        bearing = ((3 * pi / 2) - theta) % (2 * pi)
    dist = np.sqrt(u_traj ** 2 + v_traj ** 2) * d_time
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    traced_lat_rad = np.arcsin(
        np.sin(lat_rad) * np.cos(dist / r_earth)
        + np.cos(lat_rad) * np.sin(dist / r_earth) * np.cos(bearing)
    )
    traced_lon_rad = lon_rad + np.arctan2(
        np.sin(bearing) * np.sin(dist / r_earth) * np.cos(lat_rad),
        np.cos(dist / r_earth) - np.sin(lat_rad) * np.sin(lat_rad),
    )
    traced_lat = np.rad2deg(traced_lat_rad)
    traced_lon = longitude_set_meridian(np.rad2deg(traced_lon_rad))
    return traced_lat, traced_lon


def trace_forward(lat, lon, u_traj, v_traj, d_time):
    """Wrapper function for forward trajectory, here d_time is positive"""
    return trace_one_way(lat, lon, u_traj, v_traj, d_time, True)


def trace_backward(lat, lon, u_traj, v_traj, d_time):
    """Wrapper function for backward trajectory, here d_time is positive"""
    return trace_one_way(lat, lon, u_traj, v_traj, d_time, False)


def trace_two_way(lat, lon, u_traj, v_traj, d_time_forward, d_time_total):
    """Calculates both previous and next position given a point in between"""
    d_time_backward = d_time_total - d_time_forward

    forward_lat, forward_lon = trace_one_way(
        lat, lon, u_traj, v_traj, d_time_forward, True
    )
    backward_lat, backward_lon = trace_one_way(
        lat, lon, u_traj, v_traj, d_time_backward, False
    )
    return forward_lat, forward_lon, backward_lat, backward_lon


def cos_transition(absolute_input, transition_start, transition_end):
    """function that smoothly transitions from 1 to 0 using a cosine-shaped
    transition between start and end"""
    normalised_input = (absolute_input - transition_start) / (
        transition_end - transition_start
    )
    weight_factor = 1.0 * (normalised_input < 0.0) + (
        0.5 + 0.5 * np.cos(normalised_input * pi)
    ) * (1.0 - (normalised_input < 0.0) - (normalised_input > 1.0))
    return weight_factor


def weighted_velocity(ds_for_vel, trajectory_dict):
    """Weighted velociy: needs more work"""
    pres_cutoff_start = trajectory_dict["pres_cutoff_start"]
    pres_cutoff_end = trajectory_dict["pres_cutoff_end"]
    height_factor = cos_transition(
        ds_for_vel["p_f"][:, 1:, :, :].values, pres_cutoff_start, pres_cutoff_end
    )
    weights = (
        (ds_for_vel["p_h"][:, :-1, :, :].values - ds_for_vel["p_h"][:, 1:, :, :].values)
        * ds_for_vel["q"][:, 1:, :, :].values
        * height_factor
    )
    inv_weights = 1.0 / np.sum(weights)
    u_weighted = inv_weights * np.sum(ds_for_vel["u"][:, 1:, :, :].values * weights)
    v_weighted = inv_weights * np.sum(ds_for_vel["v"][:, 1:, :, :].values * weights)
    return u_weighted, v_weighted


def get_velocity_from_strategy(ds_column, trajectory_dict):
    """wrapper routine, determine velocity according to strategy"""
    if trajectory_dict["velocity_strategy"] == "lower_troposphere_humidity_weighted":
        u_traj, v_traj = weighted_velocity(ds_column, trajectory_dict)
    else:
        raise NotImplementedError("Trajectory velocity strategy not implemented")
    return u_traj, v_traj


def add_globals_attrs_to_ds(ds_to_add_to):
    """Adds global attributes to datasets"""
    global_attrs = {
        r"Conventions": r"CF-1.7",
        r"Contact": r"l.c.denby[at]leeds[dot]ac[dot again]uk s.boeing[at]leeds[dot]ac[dot again]uk",
        r"ERA5 reference": r"Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., ... & Simmons, A. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society.",
        r"Created": datetime.datetime.now().isoformat(),
        r"Created with": r"https://github.com/EUREC4A-UK/lagtraj",
    }
    for attribute in global_attrs:
        ds_to_add_to.attrs[attribute] = global_attrs[attribute]


def add_dict_to_global_attrs(ds_to_add_to, dictionary):
    """Adds global attributes to datasets"""
    for attribute in dictionary:
        ds_to_add_to.attrs[attribute] = dictionary[attribute]


def fix_units(ds_to_fix):
    """Changes units of ERA5 data to make them compatible with the cf-checker"""
    units_dict = {
        "(0 - 1)": "1",
        "m of water equivalent": "m",
        "~": "1",
    }
    for variable in ds_to_fix.variables:
        if hasattr(ds_to_fix[variable], "units"):
            these_units = ds_to_fix[variable].units
            if these_units in units_dict:
                ds_to_fix[variable].attrs["units"] = units_dict[these_units]


def stationary_trajectory(ds_traj, trajectory_dict):
    """Adds data for a target point that is directly in the time series"""
    lat_target = trajectory_dict["lat_origin"]
    lon_target = longitude_set_meridian(trajectory_dict["lon_origin"])
    ds_traj["lat_traj"][:] = lat_target
    ds_traj["lon_traj"][:] = lon_target
    ds_traj["u_traj"][:] = 0.0
    ds_traj["v_traj"][:] = 0.0
    ds_traj["processed"][:] = True


def prescribed_velocity_trajectory(ds_traj, trajectory_dict):
    """Adds data for a target point that is directly in the time series"""
    lat_target = trajectory_dict["lat_origin"]
    lon_target = longitude_set_meridian(trajectory_dict["lon_origin"])
    time_target = np.datetime64(trajectory_dict["datetime_origin"])
    u_traj = trajectory_dict["u_traj"]
    v_traj = trajectory_dict["v_traj"]
    for index in range(len(ds_traj["time"])):
        d_time = (
            (ds_traj["time"][index] - time_target) / np.timedelta64(1, "s")
        ).values
        if d_time < 0.0:
            lat_at_time, lon_at_time = trace_backward(
                lat_target, lon_target, u_traj, v_traj, -d_time
            )
        else:
            lat_at_time, lon_at_time = trace_forward(
                lat_target, lon_target, u_traj, v_traj, d_time
            )
        ds_traj["lat_traj"][index] = lat_at_time
        ds_traj["lon_traj"][index] = lon_at_time
        ds_traj["u_traj"][index] = u_traj
        ds_traj["v_traj"][index] = v_traj
        ds_traj["processed"][index] = True


def trajectory_at_target(ds_time_selection, ds_traj, trajectory_dict):
    """Adds data for a target point that is directly in the time series"""
    time_target = np.datetime64(trajectory_dict["datetime_origin"])
    lat_target = trajectory_dict["lat_origin"]
    lon_target = longitude_set_meridian(trajectory_dict["lon_origin"])
    ds_time = ds_time_selection.sel(time=[time_target])
    time_exact_index = np.argmax(ds_time_selection["time"] == time_target)
    ds_local = era5_interp_column(ds_time, lat_target, lon_target)
    add_heights_and_pressures(ds_local)
    u_traj, v_traj = get_velocity_from_strategy(ds_local, trajectory_dict)
    ds_traj["lat_traj"][time_exact_index] = lat_target
    ds_traj["lon_traj"][time_exact_index] = lon_target
    ds_traj["u_traj"][time_exact_index] = u_traj
    ds_traj["v_traj"][time_exact_index] = v_traj
    ds_traj["processed"][time_exact_index] = True


def trajectory_around_target(ds_time_selection, ds_traj, trajectory_dict):
    """Adds data around a target point that is not directly in the time series and 
    needs interpolation"""
    nr_iterations_traj = trajectory_dict["nr_iterations_traj"]
    time_target = np.datetime64(trajectory_dict["datetime_origin"])
    lat_target = trajectory_dict["lat_origin"]
    lon_target = longitude_set_meridian(trajectory_dict["lon_origin"])
    # Find relevant indices
    time_greater_index = np.argmax(ds_time_selection["time"] > time_target)
    time_smaller_index = time_greater_index - 1
    ds_interpolated = era5_interp_column(ds_time_selection, lat_target, lon_target)
    add_heights_and_pressures(ds_interpolated)
    u_guess, v_guess = get_velocity_from_strategy(ds_interpolated, trajectory_dict)
    d_time_forward = (
        (ds_time_selection["time"][time_greater_index] - time_target)
        / np.timedelta64(1, "s")
    ).values
    d_time_total = (
        (
            ds_time_selection["time"][time_greater_index]
            - ds_time_selection["time"][time_smaller_index]
        )
        / np.timedelta64(1, "s")
    ).values
    # iteratively find velocity at adjacent points in time
    for _ in range(nr_iterations_traj):
        forward_lat, forward_lon, backward_lat, backward_lon = trace_two_way(
            lat_target, lon_target, u_guess, v_guess, d_time_forward, d_time_total
        )
        ds_begin = ds_time_selection.isel(time=[time_smaller_index])
        ds_begin_column = era5_interp_column(ds_begin, backward_lat, backward_lon)
        add_heights_and_pressures(ds_begin_column)
        u_begin, v_begin = get_velocity_from_strategy(ds_begin_column, trajectory_dict)
        ds_end = ds_time_selection.isel(time=[time_greater_index])
        ds_end_column = era5_interp_column(ds_end, forward_lat, forward_lon)
        add_heights_and_pressures(ds_end_column)
        u_end, v_end = get_velocity_from_strategy(ds_end_column, trajectory_dict)
        u_guess = 0.5 * (u_begin + u_end)
        v_guess = 0.5 * (v_begin + v_end)
    ds_traj["lat_traj"][time_smaller_index] = backward_lat
    ds_traj["lon_traj"][time_smaller_index] = backward_lon
    ds_traj["u_traj"][time_smaller_index] = u_begin
    ds_traj["v_traj"][time_smaller_index] = v_begin
    ds_traj["processed"][time_smaller_index] = True
    ds_traj["lat_traj"][time_greater_index] = forward_lat
    ds_traj["lon_traj"][time_greater_index] = forward_lon
    ds_traj["u_traj"][time_greater_index] = u_end
    ds_traj["v_traj"][time_greater_index] = v_end
    ds_traj["processed"][time_greater_index] = True


def forward_trajectory(ds_time_selection, ds_traj, trajectory_dict):
    """Adds data after the last index that is already initialised"""
    # Ugly
    last_processed_index = np.argmax(
        ds_traj["processed"] * np.arange(len(ds_traj["processed"]))
    ).values
    nr_iterations_traj = trajectory_dict["nr_iterations_traj"]
    for forward_index in range(last_processed_index + 1, len(ds_traj["processed"])):
        d_time_forward = (
            (
                ds_time_selection["time"][forward_index]
                - ds_time_selection["time"][forward_index - 1]
            )
            / np.timedelta64(1, "s")
        ).values
        u_begin = ds_traj["u_traj"][forward_index - 1]
        v_begin = ds_traj["v_traj"][forward_index - 1]
        lat_begin = ds_traj["lat_traj"][forward_index - 1]
        lon_begin = ds_traj["lon_traj"][forward_index - 1]
        u_guess = u_begin
        v_guess = v_begin
        # iteratively find velocity at adjacent points in time
        for _ in range(nr_iterations_traj):
            forward_lat, forward_lon = trace_forward(
                lat_begin, lon_begin, u_guess, v_guess, d_time_forward
            )
            ds_end = ds_time_selection.isel(time=[forward_index])
            ds_end_column = era5_interp_column(ds_end, forward_lat, forward_lon)
            add_heights_and_pressures(ds_end_column)
            u_end, v_end = get_velocity_from_strategy(ds_end_column, trajectory_dict)
            u_guess = 0.5 * (u_begin + u_end)
            v_guess = 0.5 * (v_begin + v_end)
        ds_traj["lat_traj"][forward_index] = forward_lat
        ds_traj["lon_traj"][forward_index] = forward_lon
        ds_traj["u_traj"][forward_index] = u_end
        ds_traj["v_traj"][forward_index] = v_end
        ds_traj["processed"][forward_index] = True


def backward_trajectory(ds_time_selection, ds_traj, trajectory_dict):
    """Adds data before the first index that is already initialised"""
    # Ugly
    first_processed_index = np.argmax(ds_traj["processed"]).values
    nr_iterations_traj = trajectory_dict["nr_iterations_traj"]
    for backward_index in range(first_processed_index - 1, -1, -1):
        d_time_backward = (
            (
                ds_time_selection["time"][backward_index + 1]
                - ds_time_selection["time"][backward_index]
            )
            / np.timedelta64(1, "s")
        ).values
        u_end = ds_traj["u_traj"][backward_index + 1]
        v_end = ds_traj["v_traj"][backward_index + 1]
        lat_end = ds_traj["lat_traj"][backward_index + 1]
        lon_end = ds_traj["lon_traj"][backward_index + 1]
        u_guess = u_end
        v_guess = v_end
        # iteratively find velocity at adjacent points in time
        for _ in range(nr_iterations_traj):
            backward_lat, backward_lon = trace_backward(
                lat_end, lon_end, u_guess, v_guess, d_time_backward
            )
            ds_begin = ds_time_selection.isel(time=[backward_index])
            ds_begin_column = era5_interp_column(ds_begin, backward_lat, backward_lon)
            add_heights_and_pressures(ds_begin_column)
            u_begin, v_begin = get_velocity_from_strategy(
                ds_begin_column, trajectory_dict
            )
            u_guess = 0.5 * (u_begin + u_end)
            v_guess = 0.5 * (v_begin + v_end)
        ds_traj["lat_traj"][backward_index] = backward_lat
        ds_traj["lon_traj"][backward_index] = backward_lon
        ds_traj["u_traj"][backward_index] = u_begin
        ds_traj["v_traj"][backward_index] = v_begin
        ds_traj["processed"][backward_index] = True


def dummy_trajectory(mf_dataset, trajectory_dict):
    """Trajectory example: needs to be integrated into main functionality"""
    time_target = np.datetime64(trajectory_dict["datetime_origin"])
    start_date = time_target - np.timedelta64(
        trajectory_dict["backward_duration_hours"], "h"
    )
    end_date = time_target + np.timedelta64(
        trajectory_dict["forward_duration_hours"], "h"
    )
    time_start_mf = np.max(mf_dataset["time"].where(mf_dataset["time"] <= start_date))
    time_end_mf = np.min(mf_dataset["time"].where(mf_dataset["time"] >= end_date))
    ds_time_selection = mf_dataset.sel(time=slice(time_start_mf, time_end_mf))
    ds_traj = xr.Dataset(coords={"time": ds_time_selection.time})
    time_len = len(ds_time_selection["time"].values)
    ds_traj["lat_traj"] = (
        ("time"),
        np.empty((time_len)),
        {"long_name": "trajectory latitude", "units": "degrees_east"},
    )
    ds_traj["lon_traj"] = (
        ("time"),
        np.empty((time_len)),
        {"long_name": "trajectory longitude", "units": "degrees_north"},
    )
    ds_traj["u_traj"] = (
        ("time"),
        np.empty((time_len)),
        {"long_name": "trajectory U component of wind", "units": "m s**-1"},
    )
    ds_traj["v_traj"] = (
        ("time"),
        np.empty((time_len)),
        {"long_name": "trajectory V component of wind", "units": "m s**-1"},
    )
    ds_traj["processed"] = (
        ("time"),
        np.empty((time_len)),
        {"long_name": "Has array been processed", "units": "-"},
    )
    ds_traj["processed"].values[:] = False
    time_exact_match = time_target in ds_time_selection["time"]
    if trajectory_dict["velocity_strategy"] == "stationary":
        stationary_trajectory(ds_traj, trajectory_dict)
    if trajectory_dict["velocity_strategy"] == "prescribed_velocity":
        prescribed_velocity_trajectory(ds_traj, trajectory_dict)
    elif time_exact_match:
        trajectory_at_target(ds_time_selection, ds_traj, trajectory_dict)
        forward_trajectory(ds_time_selection, ds_traj, trajectory_dict)
        backward_trajectory(ds_time_selection, ds_traj, trajectory_dict)
    else:
        trajectory_around_target(ds_time_selection, ds_traj, trajectory_dict)
        forward_trajectory(ds_time_selection, ds_traj, trajectory_dict)
        backward_trajectory(ds_time_selection, ds_traj, trajectory_dict)
    if not all(ds_traj["processed"].values[:]):
        raise Exception("Trajectory issue, not all timesteps have been filled")
    ds_traj = ds_traj.drop_vars(["processed"])
    fix_units(ds_traj)
    add_globals_attrs_to_ds(ds_traj)
    add_dict_to_global_attrs(ds_traj, trajectory_dict)
    ds_traj.to_netcdf("ds_traj.nc")


def dummy_forcings(mf_dataset, forcings_dict):
    """Forcings example: needs to be integrated into main functionality"""
    ds_out = xr.Dataset()
    ds_traj = xr.open_dataset(forcings_dict["traj_file"])
    for index in range(len(ds_traj["time"])):
        this_time = ds_traj["time"][index]
        # Ugly
        mf_index = np.argmax(mf_dataset["time"] == this_time).values
        ds_time = mf_dataset.isel(time=[mf_index])
        half_averaging_width = 0.5 * forcings_dict["averaging_width"]
        lats_lons_dict = {
            "lat_min": ds_traj["lat_traj"][index].values - half_averaging_width,
            "lat_max": ds_traj["lat_traj"][index].values + half_averaging_width,
            "lon_min": longitude_set_meridian(ds_traj["lon_traj"][index].values)
            - half_averaging_width,
            "lon_max": longitude_set_meridian(ds_traj["lon_traj"][index].values)
            + half_averaging_width,
            "lat": ds_traj["lat_traj"][index].values,
            "lon": longitude_set_meridian(ds_traj["lon_traj"][index].values),
            "u_traj": ds_traj["u_traj"][index].values,
            "v_traj": ds_traj["v_traj"][index].values,
        }
        lats_lons_dict.update(forcings_dict)
        out_levels = np.arange(0, 10000.0, 40.0)
        ds_smaller = era_5_subset(ds_time, lats_lons_dict)
        add_heights_and_pressures(ds_smaller)
        add_auxiliary_variables(ds_smaller, ["theta", "rho"])
        ds_time_height = era5_on_height_levels(ds_smaller, out_levels)
        era5_add_lat_lon_meshgrid(ds_time_height)
        ds_profiles = era5_single_point(ds_time_height, lats_lons_dict)
        ds_era5_mean = era5_box_mean(ds_time_height, lats_lons_dict)
        for variable in ds_era5_mean.variables:
            if variable not in ["time", "lev"]:
                ds_profiles[variable + "_mean"] = ds_era5_mean[variable]
        ds_gradients = era5_gradients(
            ds_time_height, ["u", "v", "p_f", "theta"], lats_lons_dict
        )
        ds_time_step = xr.merge((ds_gradients, ds_profiles))
        ds_tendencies = era5_adv_tendencies(
            ds_time_step, ["u", "v", "theta"], lats_lons_dict
        )
        ds_time_step = xr.merge((ds_time_step, ds_tendencies))
        add_geowind_around_centre(ds_time_step, lats_lons_dict)
        ds_time_step.reset_coords(["latitude", "longitude"])
        ds_out = xr.combine_by_coords((ds_out, ds_time_step))
    # Add trajectory information
    ds_out = xr.combine_by_coords((ds_out, ds_traj))
    fix_units(ds_out)
    add_globals_attrs_to_ds(ds_out)
    add_dict_to_global_attrs(ds_out, forcings_dict)
    ds_out.to_netcdf("ds_along_traj.nc")


def main():
    """Dummy implementations for trajectory tool"""
    files_model_an = "output_domains/model_an_*_eurec4a_circle_eul_domain.nc"
    files_single_an = "output_domains/single_an_*_eurec4a_circle_eul_domain.nc"
    files_model_fc = "output_domains/model_fc_*_eurec4a_circle_eul_domain.nc"
    files_single_fc = "output_domains/single_fc_*_eurec4a_circle_eul_domain.nc"
    ds_model_an = xr.open_mfdataset(files_model_an, combine="by_coords")
    # z needs to be dropped to prevent duplicity, lnsp is simply redundant
    ds_model_an = ds_model_an.drop_vars(["z", "lnsp"])
    ds_single_an = xr.open_mfdataset(files_single_an, combine="by_coords")
    ds_model_fc = xr.open_mfdataset(files_model_fc, combine="by_coords")
    ds_single_fc = xr.open_mfdataset(files_single_fc, combine="by_coords")
    ds_list = [ds_model_an, ds_single_an, ds_model_fc, ds_single_fc]
    for this_ds in ds_list:
        era_5_normalise_longitude(this_ds)
    ds_merged = xr.merge(ds_list)
    dummy_trajectory_dict = {
        "lat_origin": 13.3,
        "lon_origin": -57.717,
        "datetime_origin": "2020-02-03T12:30",
        "backward_duration_hours": 3,
        "forward_duration_hours": 1,
        "nr_iterations_traj": 10,
        "velocity_strategy": "lower_troposphere_humidity_weighted",
        # "velocity_strategy": "prescribed_velocity",
        # "u_traj" : -6.0,
        # "v_traj" : -0.25,
        "pres_cutoff_start": 60000.0,
        "pres_cutoff_end": 50000.0,
    }
    dummy_trajectory(ds_merged, dummy_trajectory_dict)
    dummy_forcings_dict = {
        "gradients_strategy": "both",
        "mask": "ocean",
        "traj_file": "ds_traj.nc",
        "averaging_width": 4.0,
    }
    dummy_forcings(ds_merged, dummy_forcings_dict)


if __name__ == "__main__":
    main()
