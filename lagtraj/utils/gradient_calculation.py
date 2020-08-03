import numpy as np
import xarray as xr


from .masking import era5_mask
from .distance import lat_dist, lon_dist
from ..interpolation import boundary_gradients
from .geometry import longitude_set_meridian


def _era5_boundary_gradients(da_box):
    """ Calculate gradients from boundary values
    Using distances along latitude and longitude axes
    Weight by box size?"""
    da_filtered = da_box
    lat1_point = dictionary["lat"]
    lon1_point = longitude_set_meridian(dictionary["lon"])
    lat2_mg = da_filtered["lat_meshgrid"].values
    lon2_mg = da_filtered["lon_meshgrid"].values
    # Use the distance at the actual latitude to calculate gradients between boundaries
    x_array = lon_dist(lon1_point, lon2_mg, lat2_mg)
    y_array = lat_dist(lat1_point, lat2_mg)
    val_array = da_filtered.values
    return boundary_gradients(x_array, y_array, val_array)


def _era5_regression_gradients(da_box):
    """ Calculate gradients function using local coordinate system and
    Regression, using the normal equation"""
    da_filtered = da_box
    if "mask" in dictionary:
        mask = era5_mask(da_box, dictionary)
        da_filtered = da_filtered.where(mask)
    lat1_point = dictionary["lat"]
    lon1_point = longitude_set_meridian(dictionary["lon"])
    lat2_mg = da_filtered["lat_meshgrid"].values
    lon2_mg = da_filtered["lon_meshgrid"].values
    # Use the center latitude for projection onto plane
    x_array = lon_dist(lon1_point, lon2_mg, lat1_point)
    y_array = lat_dist(lat1_point, lat2_mg)
    x_flat = x_array.flatten()
    y_flat = y_array.flatten()
    ones_flat = np.ones(np.shape(x_flat))
    len_temp = np.shape(da_box)[0]
    len_levels = np.shape(da_box)[1]
    x_gradient_array = np.empty((len_temp, len_levels))
    y_gradient_array = np.empty((len_temp, len_levels))
    for this_time in range(len_temp):
        for this_level in range(len_levels):
            data_flat = da_filtered[
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


def calc_horizontal_gradients(da_field, method="regression"):
    """
    Compute horizontal gradients in `da_field` using `method`. Any nan-values
    in `da_field` will be ignored in the gradient calculation
    """

    if method == "regression":
        x_gradient_array, y_gradient_array = _era5_regression_gradients(
            da_field
        )
    elif method == "boundary":
        x_gradient_array, y_gradient_array = _era5_boundary_gradients(
            da_field
        )
    else:
        raise NotImplementedError("Gradients strategy not implemented")

    v = da_field.name
    da_dphidx = xr.DataArray(
        x_gradient_array,
        name=f"d{v}dx",
        dims=("time", "lev"),
        attrs=dict(
            long_name=f"f{da_field.long_name} x-gradient",
            units=f"{da_field.units} m**-1",
        )
    )
    da_dphidy = xr.DataArray(
        y_gradient_array,
        name=f"d{v}dy",
        dims=("time", "lev"),
        attrs=dict(
            long_name=f"f{da_field.long_name} y-gradient",
            units=f"{da_field.units} m**-1",
        )
    )

    return da_dphidx, da_dphidy
