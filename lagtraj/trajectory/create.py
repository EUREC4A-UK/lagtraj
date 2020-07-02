import datetime
from pathlib import Path

import xarray as xr
import numpy as np

from .. import DEFAULT_ROOT_DATA_PATH
from .load import load_definition
from . import build_data_path
from ..domain.load import load_data as load_domain_data
from ..utils import optional_debugging

""" Routines for creating a trajectory

# - Implement different strategies (single height, weighted over heights, in
#   future possibly hysplit)
- Add metadata to NetCDF output
- Improve linear trajectory to work with haversine functions and actual velocities
- Relax assumption of hourly data?"""


def create_trajectory(origin, trajectory_type, da_times, **kwargs):
    trajectory_fn = None
    if trajectory_type == "eulerian":
        trajectory_fn = create_eulerian_trajectory
    elif trajectory_type == "linear":
        trajectory_fn = create_linear_trajectory
    elif trajectory_type == "single_level":
        trajectory_fn = create_single_level_trajectory
    elif trajectory_type == "weighted":
        trajectory_fn = create_weighted_trajectory
    else:
        raise Exception("Trajectory_type not found")

    return trajectory_fn(origin=origin, da_times=da_times, **kwargs)


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("trajectory")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument("--debug", default=False, action="store_true")
    args = argparser.parse_args()

    cli(
        data_path=args.data_path, trajectory_name=args.trajectory_name, debug=args.debug
    )


def cli(data_path, trajectory_name, debug):
    traj_definition = load_definition(root_data_path=data_path, name=trajectory_name)

    if traj_definition.timestep == "domain_data":
        da_times = _get_times_from_domain(
            trajectory_definition=traj_definition, root_data_path=data_path
        )
    elif type(traj_definition.timestep) == datetime.timedelta:
        da_times = _build_times_dataarray(
            origin=traj_definition.origin,
            duration=traj_definition.duration,
            dt=traj_definition.timestep,
        )
    else:
        raise NotImplementedError(traj_definition.timestep)

    with optional_debugging(debug):
        ds_trajectory = create_trajectory(
            origin=traj_definition.origin,
            trajectory_type=traj_definition.type,
            da_times=da_times,
        )

    trajectory_data_path = build_data_path(
        root_data_path=data_path, trajectory_name=traj_definition.name
    )

    ds_trajectory.to_netcdf(trajectory_data_path)
    print("Saved trajectory to `{}`".format(trajectory_data_path))


def _get_times_from_domain(trajectory_definition, root_data_path):
    ds_domain = load_domain_data(
        root_data_path=root_data_path, name=trajectory_definition.domain
    )
    t0 = trajectory_definition.origin.datetime
    t_min = t0 - trajectory_definition.duration.backward
    t_max = t0 + trajectory_definition.duration.forward
    da_times = ds_domain.sel(time=slice(t_min, t_max)).time
    if da_times.count() == 0:
        raise Exception(
            "You selected to use the domain data for timesteps"
            " in the trajectory, but in the time interval selected"
            " for the trajectory ({}, {}) there are is no domain data"
            " (time range: {} to {})".format(
                t_min,
                t_max,
                ds_domain.time.min().dt.strftime("%Y-%m-%d %H:%M").item(),
                ds_domain.time.max().dt.strftime("%Y-%m-%d %H:%M").item(),
            )
        )
    return da_times


def _build_times_dataarray(origin, duration, dt):
    t0 = origin.datetime
    t_min = origin.datetime + duration.backward
    t_max = origin.datetime + duration.forward

    times = [t0]
    while times[0] > t_min:
        t_new = times[-1] + dt
        times.insert(0, t_new)
    while times[-1] < t_max:
        t_new = times[-1] + dt
        times.append(t_new)

    return xr.DataArray(times, name="time", dims=("time"))


def create_eulerian_trajectory(origin, da_times, ds_domain=None):
    ds = xr.Dataset(coords=dict(time=da_times))

    lat0 = origin.lat
    lon0 = origin.lon
    ds["lat"] = (
        ("time",),
        lat0 * np.ones(len(ds.time)),
        {"long_name": "latitude", "units": "degrees_east"},
    )
    ds["lon"] = (
        ("time",),
        lon0 * np.ones(len(ds.time)),
        {"long_name": "longitude", "units": "degrees_north"},
    )

    return ds


def create_linear_trajectory(duration, origin, ds_domain=None):
    raise NotImplementedError

    # times = pd.date_range(
    #     np.datetime64(trajectory_params["datetime_end"])
    #     - np.timedelta64(trajectory_params["duration_hours"], "h"),
    #     periods=trajectory_params["duration_hours"] + 1,
    #     freq="h",
    # )
    # nr_hours = len(times)
    # lat_end = trajectory_params["lat_end"]
    # lon_end = trajectory_params["lon_end"]
    # dlat_dt = trajectory_params["dlat_dt"]
    # dlon_dt = trajectory_params["dlon_dt"]
    # lats = lat_end - 3600.0 * dlat_dt * np.arange(nr_hours - 1, -1, -1)
    # lons = lon_end - 3600.0 * dlon_dt * np.arange(nr_hours - 1, -1, -1)
    # return trajectory_to_xarray(times, lats, lons)


def stationary_trajectory(ds_traj, trajectory_dict):
    """Adds data for a stationary origin point"""
    lat_origin = trajectory_dict["lat_origin"]
    lon_origin = longitude_set_meridian(trajectory_dict["lon_origin"])
    ds_traj["lat_traj"][:] = lat_origin
    ds_traj["lon_traj"][:] = lon_origin
    ds_traj["u_traj"][:] = 0.0
    ds_traj["v_traj"][:] = 0.0
    ds_traj["processed"][:] = True


def prescribed_velocity_trajectory(ds_traj, trajectory_dict):
    """Adds data for origin point using constant velocity"""
    lat_origin = trajectory_dict["lat_origin"]
    lon_origin = longitude_set_meridian(trajectory_dict["lon_origin"])
    time_origin = np.datetime64(trajectory_dict["datetime_origin"])
    u_traj = trajectory_dict["u_traj"]
    v_traj = trajectory_dict["v_traj"]
    for index in range(len(ds_traj["time"])):
        d_time = (
            (ds_traj["time"][index] - time_origin) / np.timedelta64(1, "s")
        ).values
        if d_time < 0.0:
            lat_at_time, lon_at_time = trace_backward(
                lat_origin, lon_origin, u_traj, v_traj, -d_time
            )
        else:
            lat_at_time, lon_at_time = trace_forward(
                lat_origin, lon_origin, u_traj, v_traj, d_time
            )
        ds_traj["lat_traj"][index] = lat_at_time
        ds_traj["lon_traj"][index] = lon_at_time
        ds_traj["u_traj"][index] = u_traj
        ds_traj["v_traj"][index] = v_traj
        ds_traj["processed"][index] = True



def dummy_trajectory(mf_dataset, trajectory_dict):
    """Trajectory example: needs to be integrated into main functionality"""
    time_origin = np.datetime64(trajectory_dict["datetime_origin"])
    start_date = time_origin - np.timedelta64(
        trajectory_dict["backward_duration_hours"], "h"
    )
    end_date = time_origin + np.timedelta64(
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
    time_exact_match = time_origin in ds_time_selection["time"]

    trajectory_at_origin(ds_time_selection, ds_traj, trajectory_dict)
    forward_trajectory(ds_time_selection, ds_traj, trajectory_dict)
    backward_trajectory(ds_time_selection, ds_traj, trajectory_dict)

    if not all(ds_traj["processed"].values[:]):
        raise Exception("Trajectory issue, not all timesteps have been filled")
    ds_traj = ds_traj.drop_vars(["processed"])
    fix_units(ds_traj)
    add_globals_attrs_to_ds(ds_traj)
    add_dict_to_global_attrs(ds_traj, trajectory_dict)
    ds_traj.to_netcdf("ds_traj.nc")


def create_single_level_trajectory(duration, origin):
    raise NotImplementedError


def create_weighted_trajectory(duration, origin):
    raise NotImplementedError


if __name__ == "__main__":
    main()
