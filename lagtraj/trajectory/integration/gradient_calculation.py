import numpy as np
import xarray as xr


from .masking import era5_mask
from .distance import lat_dist, lon_dist
from ..interpolation import boundary_gradients
from .geometry import longitude_set_meridian


def _era5_boundary_gradients(ds_box, variable, dictionary):
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
    # Use the distance at the actual latitude to calculate gradients between boundaries
    x_array = lon_dist(lon1_point, lon2_mg, lat2_mg)
    y_array = lat_dist(lat1_point, lat2_mg)
    val_array = ds_filtered[variable].values
    return boundary_gradients(x_array, y_array, val_array)


def _era5_regression_gradients(ds_box, variable, dictionary):
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
    # Use the center latitude for projection onto plane
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
            x_gradient_array, y_gradient_array = _era5_regression_gradients(
                ds_field, variable, dictionary
            )
        elif dictionary["gradients_strategy"] == "boundary":
            x_gradient_array, y_gradient_array = _era5_boundary_gradients(
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
            x_gradient_array, y_gradient_array = _era5_boundary_gradients(
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
