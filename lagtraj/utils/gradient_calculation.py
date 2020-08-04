import numpy as np
import xarray as xr


from .. import njit
from .distance import lat_dist, lon_dist


@njit
def _boundary_gradients(x_array, y_array, val_array):
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
                vals_at_lat = val_array[this_time, this_level, this_lat, :].flatten()
                x_at_lat = x_array[this_lat, :].flatten()
                vals_filtered = vals_at_lat[~np.isnan(vals_at_lat)]
                x_filtered = x_at_lat[~np.isnan(vals_at_lat)]
                dvals = vals_filtered[-1] - vals_filtered[0]
                dval_dx[this_lat] = dvals / (x_filtered[-1] - x_filtered[0])
            for this_lon in range(len_lons):
                vals_at_lon = val_array[this_time, this_level, :, this_lon].flatten()
                y_at_lat = y_array[:, this_lon].flatten()
                vals_filtered = vals_at_lon[~np.isnan(vals_at_lon)]
                y_filtered = y_at_lat[~np.isnan(vals_at_lon)]
                dvals = vals_filtered[-1] - vals_filtered[0]
                dval_dy[this_lon] = dvals / (y_filtered[-1] - y_filtered[0])
            x_gradient_array[this_time, this_level] = np.mean(
                dval_dx[~np.isnan(dval_dx)]
            )
            y_gradient_array[this_time, this_level] = np.mean(
                dval_dy[~np.isnan(dval_dy)]
            )
    return x_gradient_array, y_gradient_array


def _era5_boundary_gradients(da_box, ds_ref_pt):
    """ Calculate gradients from boundary values
    Using distances along latitude and longitude axes
    Weight by box size?"""
    lat1_point, lon1_point = ds_ref_pt.lat.values, ds_ref_pt.lon.values
    lon2_mg, lat2_mg = np.meshgrid(da_box.lon.values, da_box.lat.values)
    # Use the distance at the actual latitude to calculate gradients between boundaries
    x_array = lon_dist(lon1_point, lon2_mg, lat2_mg)
    y_array = lat_dist(lat1_point, lat2_mg)
    val_array = da_box.values
    return _boundary_gradients(x_array, y_array, val_array)


def _era5_regression_gradients(da_box, ds_ref_pt):
    """ Calculate gradients function using local coordinate system and
    Regression, using the normal equation"""
    lat1_point, lon1_point = ds_ref_pt.lat.values, ds_ref_pt.lon.values
    lon2_mg, lat2_mg = np.meshgrid(da_box.lon.values, da_box.lat.values)
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
            data_flat = da_box[this_time, this_level, :, :].values.flatten()
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


def calc_horizontal_gradients(da_field, ds_ref_pt, method="regression"):
    """
    Compute horizontal gradients in `da_field` using `method`. Any nan-values
    in `da_field` will be ignored in the gradient calculation
    """

    if method == "regression":
        x_gradient_array, y_gradient_array = _era5_regression_gradients(
            da_field, ds_ref_pt=ds_ref_pt
        )
    elif method == "boundary":
        x_gradient_array, y_gradient_array = _era5_boundary_gradients(
            da_field, ds_ref_pt=ds_ref_pt
        )
    else:
        raise NotImplementedError(f"Gradient method `{method}` not implemented")

    v = da_field.name
    da_dphidx = xr.DataArray(
        x_gradient_array,
        name=f"d{v}dx",
        dims=("time", "lev"),
        attrs=dict(
            long_name=f"{da_field.long_name} x-gradient",
            units=f"{da_field.units} m**-1",
        ),
    )
    da_dphidy = xr.DataArray(
        y_gradient_array,
        name=f"d{v}dy",
        dims=("time", "lev"),
        attrs=dict(
            long_name=f"{da_field.long_name} y-gradient",
            units=f"{da_field.units} m**-1",
        ),
    )

    return da_dphidx, da_dphidy
