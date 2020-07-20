import numpy as np
from scipy.constants import pi

from ..utils.constants import r_earth
from ..utils import geometry

from .integration.velocity_estimation import estimate_horizontal_velocities

# optimisation constants
two_pi = 2.0 * pi
half_pi = 0.5 * pi


def extrapolate_posn_with_fixed_velocity(lat, lon, u_vel, v_vel, dt):
    """Extrapolates a position from point (lat, lon) given velocities
    (u_vel, v_vel) (in m/s) and a duration dt (in seconds)"""
    # ensure that we're just working with floats in case xarray DataArray
    # values (with coordinates) are passed in
    lat = np.float64(lat)
    lon = np.float64(lon)
    u_vel = np.float64(u_vel)
    v_vel = np.float64(v_vel)
    dt = np.float64(dt)

    if dt < 0.0:
        raise Exception("Expecting positive dt in extrapolation")
    theta = np.arctan2(v_vel, u_vel) % two_pi
    bearing = (half_pi - theta) % two_pi
    dist = np.sqrt(u_vel * u_vel + v_vel * v_vel) * dt
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
    traced_lon = geometry.longitude_set_meridian(np.rad2deg(traced_lon_rad))
    return traced_lat, traced_lon


def _extract_column_at_time(ds_domain, t, lat, lon, method):
    if method == "nearest":
        return ds_domain.sel(time=t, lat=lat, lon=lon, method="nearest")
    else:
        interp_kwargs = dict(bounds_error=True)
        return ds_domain.interp(time=t, lat=lat, lon=lon, kwargs=interp_kwargs)


def extrapolate_using_domain_data(
    lat,
    lon,
    t0,
    dt,
    ds_domain,
    velocity_method,
    velocity_method_kwargs={},
    num_velocity_integrations=10,
):
    """Using the fields in `ds_domain` apply `velocity_method` to iteratively
    estimate a velocity field and extrapolate using these velocities until the
    estimated velocities at the start and end point converge

    TODO: check for convergence"""
    time_space_interpolation = velocity_method_kwargs.pop(
        "time_space_interpolation", "linear"
    )
    ds_column_interpolated = _extract_column_at_time(
        t=t0, lat=lat, lon=lon, method=time_space_interpolation, ds_domain=ds_domain
    )

    lat_start, lon_start = lat, lon
    u_start, v_start = estimate_horizontal_velocities(
        ds_column=ds_column_interpolated,
        method=velocity_method,
        **velocity_method_kwargs
    )
    # the extrapolation function requires positive time increment for the
    # integration, so we simply reverse the time and velocities
    if dt < 0:
        s = -1.0
    else:
        s = 1.0

    u_guess, v_guess = u_start, v_start
    for n in range(num_velocity_integrations):
        if np.any(np.isnan([u_guess, v_guess])):
            raise Exception
        lat_end, lon_end = extrapolate_posn_with_fixed_velocity(
            lat=lat_start,
            lon=lon_start,
            u_vel=s * u_guess,
            v_vel=s * v_guess,
            dt=s * dt,
        )
        ds_column_interpolated = _extract_column_at_time(
            t=t0,
            lat=lat_end,
            lon=lon_end,
            method=time_space_interpolation,
            ds_domain=ds_domain,
        )
        u_end, v_end = estimate_horizontal_velocities(
            ds_column=ds_column_interpolated,
            method=velocity_method,
            **velocity_method_kwargs
        )
        u_guess = 0.5 * (u_start + u_end)
        v_guess = 0.5 * (v_start + v_end)

    return lat_end, lon_end
