import numpy as np
import xarray as xr

from .. import njit
from .geometry import lat_dist, lon_dist


@njit
def _boundary_gradients(x_array, y_array, val_array):
    """Numba function to calculate gradients, dismissing filtered points"""
    len_temp = np.shape(val_array)[0]
    len_levels = np.shape(val_array)[1]
    len_lats = np.shape(val_array)[2]
    len_lons = np.shape(val_array)[3]
    x_gradient_array = np.empty((len_temp, len_levels))
    y_gradient_array = np.empty((len_temp, len_levels))
    for this_time in range(len_temp):
        for this_level in range(len_levels):
            dxval_tot = 0.0
            dyval_tot = 0.0
            dx_tot = 0.0
            dy_tot = 0.0
            # This calculates the x-gradient as a weighted average over latitude
            # The way the averaging is done now is such that the weight of a particular
            # latitude/longitude is proportional to the length of the segment over
            # which a gradient is calculated at that latitude/longitude.
            # This length varies, and may even be zero in extreme cases,
            # due to the use of a mask (and also a bit due to the lat-lon grid,
            # but that would mostly be notable close to the poles). This changes
            # ensures there is no inappropriately high weighting given when a
            # gradient is calculated over a short (or zero) distance.
            # The ~np.isnan operation is used to filter out masked values
            # The first and last filtered values are used
            for this_lat in range(len_lats):
                vals_at_lat = val_array[this_time, this_level, this_lat, :].flatten()
                x_at_lat = x_array[this_lat, :].flatten()
                vals_filtered = vals_at_lat[~np.isnan(vals_at_lat)]
                if len(vals_filtered) > 1:
                    x_filtered = x_at_lat[~np.isnan(vals_at_lat)]
                    dxval_tot = dxval_tot + vals_filtered[-1] - vals_filtered[0]
                    dx_tot = dx_tot + x_filtered[-1] - x_filtered[0]
            # This similarly calculates the y-gradient weighted average
            for this_lon in range(len_lons):
                vals_at_lon = val_array[this_time, this_level, :, this_lon].flatten()
                y_at_lat = y_array[:, this_lon].flatten()
                vals_filtered = vals_at_lon[~np.isnan(vals_at_lon)]
                if len(vals_filtered) > 1:
                    y_filtered = y_at_lat[~np.isnan(vals_at_lon)]
                    dyval_tot = dyval_tot + vals_filtered[-1] - vals_filtered[0]
                    dy_tot = dy_tot + y_filtered[-1] - y_filtered[0]
            # Average these gradients (not weighted at this point, but filtering out all nan values due to e.g. division by zero!)
            if abs(dx_tot) > 1e-4:
                x_gradient_array[this_time, this_level] = dxval_tot / dx_tot
            else:
                x_gradient_array[this_time, this_level] = np.nan
            if abs(dy_tot) > 1e-4:
                y_gradient_array[this_time, this_level] = dyval_tot / dy_tot
            else:
                y_gradient_array[this_time, this_level] = np.nan
    return x_gradient_array, y_gradient_array


@njit
def _regression_gradients(x_array, y_array, val_array):
    """Numba function for regression gradients"""
    len_temp = np.shape(val_array)[0]
    len_levels = np.shape(val_array)[1]
    x_gradient_array = np.empty((len_temp, len_levels))
    y_gradient_array = np.empty((len_temp, len_levels))
    x_flat = x_array.flatten()
    y_flat = y_array.flatten()
    ones_flat = np.ones(np.shape(x_flat))
    for this_time in range(len_temp):
        for this_level in range(len_levels):
            # For each level and time, put all valid (non-nan) x and y locations
            # as well as the corresponding data, in 1D arrays
            data_flat = val_array[this_time, this_level, :, :].flatten()
            data_flat_filter = np.expand_dims(data_flat[~np.isnan(data_flat)], axis=1)
            x_flat_filter = np.expand_dims(x_flat[~np.isnan(data_flat)], axis=1)
            y_flat_filter = np.expand_dims(y_flat[~np.isnan(data_flat)], axis=1)
            if (np.nanmin(x_flat_filter) < np.nanmax(x_flat_filter)) and (
                np.nanmin(y_flat_filter) < np.nanmax(y_flat_filter)
            ):
                ones_flat_filter = np.expand_dims(
                    ones_flat[~np.isnan(data_flat)], axis=1
                )
                oxy_mat = np.hstack((ones_flat_filter, x_flat_filter, y_flat_filter))
                # Use the normal method to find the best fit of a plane through the data
                # At each individual model level
                theta = np.dot(
                    np.dot(
                        np.linalg.pinv(np.dot(oxy_mat.transpose(), oxy_mat)),
                        oxy_mat.transpose(),
                    ),
                    data_flat_filter,
                )
                x_gradient_array[this_time, this_level] = theta[1][0]
                y_gradient_array[this_time, this_level] = theta[2][0]
            else:
                x_gradient_array[this_time, this_level] = np.nan
                y_gradient_array[this_time, this_level] = np.nan
    return x_gradient_array, y_gradient_array


def _era5_boundary_gradients(da_box, ds_ref_pt):
    """Calculate gradients from boundary values
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
    """Calculate gradients function using local coordinate system and
    Regression, using the normal equation"""
    lat1_point, lon1_point = ds_ref_pt.lat.values, ds_ref_pt.lon.values
    lon2_mg, lat2_mg = np.meshgrid(da_box.lon.values, da_box.lat.values)
    # Use the center latitude for projection onto plane
    x_array = lon_dist(lon1_point, lon2_mg, lat1_point)
    y_array = lat_dist(lat1_point, lat2_mg)
    val_array = da_box.values
    return _regression_gradients(x_array, y_array, val_array)


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
        dims=("time", "level"),
        coords=dict(time=da_field.time, level=da_field.level),
        attrs=dict(
            long_name=f"{da_field.long_name} x-gradient",
            units=f"{da_field.units} m**-1",
        ),
    )
    da_dphidy = xr.DataArray(
        y_gradient_array,
        name=f"d{v}dy",
        dims=("time", "level"),
        coords=dict(time=da_field.time, level=da_field.level),
        attrs=dict(
            long_name=f"{da_field.long_name} y-gradient",
            units=f"{da_field.units} m**-1",
        ),
    )

    return da_dphidx, da_dphidy
