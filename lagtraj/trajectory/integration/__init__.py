import numpy as np
from scipy.constants import pi

from .velocity_estimation import get_velocity_from_strategy
from .constants import *  # shouldn't to this, make explicit what is used


def trace_one_way(lat, lon, u_traj, v_traj, d_time, lforward=True):
    """Calculates previous position given lat,lon,u_traj,v_traj, and d_time
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


def trajectory_at_origin(ds_time_selection, ds_traj, trajectory_dict):
    """Adds data for an origin point that is directly in the time series"""
    time_origin = np.datetime64(trajectory_dict["datetime_origin"])
    lat_origin = trajectory_dict["lat_origin"]
    lon_origin = longitude_set_meridian(trajectory_dict["lon_origin"])
    ds_time = ds_time_selection.sel(time=[time_origin])
    time_exact_index = np.argmax(ds_time_selection["time"] == time_origin)
    ds_local = era5_interp_column(ds_time, lat_origin, lon_origin)
    add_heights_and_pressures(ds_local)
    u_traj, v_traj = get_velocity_from_strategy(ds_local, trajectory_dict)
    ds_traj["lat_traj"][time_exact_index] = lat_origin
    ds_traj["lon_traj"][time_exact_index] = lon_origin
    ds_traj["u_traj"][time_exact_index] = u_traj
    ds_traj["v_traj"][time_exact_index] = v_traj
    ds_traj["processed"][time_exact_index] = True


def trajectory_around_origin(ds_time_selection, ds_traj, trajectory_dict):
    """Adds data around origin point that is not directly in the time series and 
    needs interpolation"""
    nr_iterations_traj = trajectory_dict["nr_iterations_traj"]
    time_origin = np.datetime64(trajectory_dict["datetime_origin"])
    lat_origin = trajectory_dict["lat_origin"]
    lon_origin = longitude_set_meridian(trajectory_dict["lon_origin"])
    # Find relevant indices
    time_greater_index = np.argmax(ds_time_selection["time"] > time_origin)
    time_smaller_index = time_greater_index - 1
    ds_interpolated = era5_time_interp_column(
        ds_time_selection, time_origin, lat_origin, lon_origin
    )
    add_heights_and_pressures(ds_interpolated)
    u_guess, v_guess = get_velocity_from_strategy(ds_interpolated, trajectory_dict)
    d_time_forward = (
        (ds_time_selection["time"][time_greater_index] - time_origin)
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
            lat_origin, lon_origin, u_guess, v_guess, d_time_forward, d_time_total
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


def integrate_trajectory(mf_list, trajectory_dict):
    """Trajectory example: needs to be integrated into main functionality"""
    time_origin = np.datetime64(trajectory_dict["datetime_origin"])
    start_date = time_origin - np.timedelta64(
        trajectory_dict["backward_duration_hours"], "h"
    )
    end_date = time_origin + np.timedelta64(
        trajectory_dict["forward_duration_hours"], "h"
    )
    mf_extract_time = []
    # Only merge time arrays for in order to fill trajectory
    # Note: may want to add some checks here so time is the same
    # between files
    for mf_list_element in mf_list:
        time_start_mf = np.max(
            mf_list_element["time"].where(mf_list_element["time"] <= start_date)
        )
        time_end_mf = np.min(
            mf_list_element["time"].where(mf_list_element["time"] >= end_date)
        )
        mf_lists_ds = xr.Dataset()
        mf_lists_ds["time"] = mf_list_element["time"].sel(
            time=slice(time_start_mf, time_end_mf)
        )
        mf_extract_time.append(mf_lists_ds)
    ds_time_selection = xr.merge(mf_extract_time)
    ds_time_selection.load()
    del mf_extract_time
    ds_traj = xr.Dataset(coords={"time": ds_time_selection.time})
    ds_time_selection.close()
    time_len = len(ds_traj["time"].values)
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
    time_exact_match = time_origin in ds_traj["time"]
    vars_for_traj = ["u", "v", "sp", "z", "t", "q", "lsm"]
    if trajectory_dict["velocity_strategy"] == "stationary":
        stationary_trajectory(ds_traj, trajectory_dict)
    if trajectory_dict["velocity_strategy"] == "prescribed_velocity":
        prescribed_velocity_trajectory(ds_traj, trajectory_dict)
    elif time_exact_match:
        trajectory_at_origin(mf_list, vars_for_traj, ds_traj, trajectory_dict)
        forward_trajectory(mf_list, vars_for_traj, ds_traj, trajectory_dict)
        backward_trajectory(mf_list, vars_for_traj, ds_traj, trajectory_dict)
    else:
        trajectory_around_origin(mf_list, vars_for_traj, ds_traj, trajectory_dict)
        forward_trajectory(mf_list, vars_for_traj, ds_traj, trajectory_dict)
        backward_trajectory(mf_list, vars_for_traj, ds_traj, trajectory_dict)
    if not all(ds_traj["processed"].values[:]):
        raise Exception("Trajectory issue, not all timesteps have been filled")
    ds_traj = ds_traj.drop_vars(["processed"])
    fix_units(ds_traj)
    add_globals_attrs_to_ds(ds_traj)
    add_dict_to_global_attrs(ds_traj, trajectory_dict)
    ds_traj.to_netcdf("ds_traj.nc")
