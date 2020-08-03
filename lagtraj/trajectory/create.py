import datetime
from pathlib import Path

import xarray as xr
import numpy as np
import tqdm

from .. import DEFAULT_ROOT_DATA_PATH
from .load import load_definition
from . import build_data_path, extrapolation
from ..domain.load import load_data as load_domain_data
from ..domain.download import download_complete
from ..utils import optional_debugging

""" Routines for creating a trajectory

# - Implement different strategies (single height, weighted over heights, in
#   future possibly hysplit)
- Add metadata to NetCDF output
- Improve linear trajectory to work with haversine functions and actual velocities
- Relax assumption of hourly data?
# fix_units(ds_traj)
# add_globals_attrs_to_ds(ds_traj)
# add_dict_to_global_attrs(ds_traj, trajectory_dict)
# ds_traj.to_netcdf("ds_traj.nc")
"""


def create_trajectory(origin, trajectory_type, da_times, **kwargs):
    if trajectory_type == "eulerian":
        return create_eulerian_trajectory(origin=origin, da_times=da_times)
    elif trajectory_type == "linear":
        if "U" not in kwargs:
            raise Exception(
                "To use the `linear` trajectory integration you"
                " must provide a velocity `U`"
            )
        return create_linear_trajectory(origin=origin, da_times=da_times, **kwargs)
    elif trajectory_type == "lagrangian":
        if "ds_domain" not in kwargs:
            raise Exception(
                "To integrate a trajectory using velocities from model data"
                " you must provide the `ds_domain` kwargs"
            )
        if "velocity_method" not in kwargs:
            raise Exception(
                "To integrate a trajectory using velocities from model data"
                " you must select a `velocity_method`"
            )
        return create_lagrangian_trajectory(origin=origin, da_times=da_times, **kwargs)
    else:
        raise NotImplementedError(f"`{trajectory_type}` trajectory type not available")


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("trajectory")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument("--debug", default=False, action="store_true")
    args = argparser.parse_args()

    with optional_debugging(args.debug):
        cli(data_path=args.data_path, trajectory_name=args.trajectory)


def cli(data_path, trajectory_name):
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

    kwargs = dict(
        origin=traj_definition.origin,
        trajectory_type=traj_definition.type,
        da_times=da_times,
        # TODO: work out how to pass these once we know which we will be using
        **traj_definition.extra_kwargs,
    )

    if traj_definition.type == "lagrangian":
        ds_domain = load_domain_data(
            root_data_path=data_path, name=traj_definition.domain
        )
        kwargs["ds_domain"] = ds_domain

    ds_trajectory = create_trajectory(**kwargs)

    trajectory_data_path = build_data_path(
        root_data_path=data_path, trajectory_name=traj_definition.name
    )

    ds_trajectory.attrs["name"] = trajectory_name
    ds_trajectory.attrs["domain_name"] = traj_definition.domain
    ds_trajectory.attrs["trajectory_type"] = traj_definition.type
    for k, v in traj_definition.extra_kwargs.items():
        if type(v) == dict:
            for k_, v_ in v.items():
                ds_trajectory.attrs[f"{k}_{k_}"] = v_
        else:
            ds_trajectory.attrs[k] = str(v)

    ds_trajectory.to_netcdf(trajectory_data_path)
    print("Saved trajectory to `{}`".format(trajectory_data_path))


def _get_times_from_domain(trajectory_definition, root_data_path):
    if not download_complete(
        root_data_path=root_data_path, domain_name=trajectory_definition.domain
    ):
        raise Exception(
            "Some of the data for the selected domain"
            f" ({trajectory_definition.domain}) hasn't been"
            " downloaded yet"
        )
    ds_domain = load_domain_data(
        root_data_path=root_data_path, name=trajectory_definition.domain
    )
    t0 = trajectory_definition.origin.datetime
    t_min = t0 - trajectory_definition.duration.backward
    t_max = t0 + trajectory_definition.duration.forward
    da_times = ds_domain.time.sel(time=slice(t_min, t_max))
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
    t_min = origin.datetime - duration.backward
    t_max = origin.datetime + duration.forward

    times = [t0]
    while times[0] > t_min:
        t_new = times[0] - dt
        times.insert(0, t_new)
    while times[-1] < t_max:
        t_new = times[-1] + dt
        times.append(t_new)

    return xr.DataArray(times, name="time", dims=("time"), coords=dict(time=times))


def create_eulerian_trajectory(origin, da_times):
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
    ds["u_traj"] = (
        ("time",),
        np.zeros(len(ds.time)),
        {"long_name": "zonal velocity", "units": "m/s"},
    )
    ds["v_traj"] = (
        ("time",),
        np.zeros(len(ds.time)),
        {"long_name": "meridional velocity", "units": "m/s"},
    )

    return ds


def create_linear_trajectory(origin, da_times, U):
    """Create linear trajectory from origin point using constant velocity"""

    def extrapolation_func(lat, lon, t0, dt):
        if dt > 0:
            s = 1.0
        else:
            s = -1.0

        lat_new, lon_new = extrapolation.extrapolate_posn_with_fixed_velocity(
            lat=lat, lon=lon, u_vel=s * U[0], v_vel=s * U[1], dt=s * dt,
        )
        u_start_and_end = (U[0], U[0])
        v_start_and_end = (U[1], U[1])
        return lat_new, lon_new, u_start_and_end, v_start_and_end

    return _create_extrapolated_trajectory(
        origin=origin, da_times=da_times, extrapolation_func=extrapolation_func
    )


def create_lagrangian_trajectory(
    origin, da_times, ds_domain, velocity_method, velocity_method_kwargs={}
):
    """Create trajectory from origin point using extracting the velocity field
    from domain data"""

    def extrapolation_func(lat, lon, t0, dt):
        return extrapolation.extrapolate_using_domain_data(
            lat=lat,
            lon=lon,
            dt=dt,
            ds_domain=ds_domain,
            t0=t0,
            velocity_method=velocity_method,
            velocity_method_kwargs=velocity_method_kwargs,
        )

    return _create_extrapolated_trajectory(
        origin=origin, da_times=da_times, extrapolation_func=extrapolation_func
    )


def _create_extrapolated_trajectory(origin, da_times, extrapolation_func):
    ds_start_posn = xr.Dataset(coords=dict(time=origin.datetime))
    ds_start_posn["lat"] = origin.lat
    ds_start_posn["lon"] = origin.lon

    da_times_backward = da_times.sel(time=slice(None, origin.datetime))
    da_times_forward = da_times.sel(time=slice(origin.datetime, None))

    # xarray doesn't have a `total_seconds` accessor for datetime objects yet
    def _calculate_seconds(da):
        return da.dt.seconds + da.dt.days * 24 * 60 * 60

    points = [ds_start_posn]

    for dir in ["backward", "forward"]:
        if dir == "backward":
            da_integrate_times = da_times_backward.values[::-1]
        elif dir == "forward":
            da_integrate_times = da_times_forward.values
        else:
            raise Exception

        for t in tqdm.tqdm(da_integrate_times):
            ds_prev_posn = points[-1]
            dt = _calculate_seconds(t - ds_prev_posn.time)
            if int(dt) == 0:
                continue
            lat, lon, u_traj, v_traj = extrapolation_func(
                lat=points[-1].lat, lon=points[-1].lon, dt=dt, t0=ds_prev_posn.time
            )
            # u_traj and v_traj are tuples containing the start and end
            # velocities for the trajectory
            u_start, u_end = u_traj
            v_start, v_end = v_traj
            # for the very first point we also want to record the starting
            # velocity which was calculated by the interpolation routine
            if len(points) == 1:
                points[0]["u_traj"] = u_start
                points[0]["v_traj"] = u_start

            ds_next_posn = xr.Dataset(coords=dict(time=t))
            ds_next_posn["lat"] = lat
            ds_next_posn["lon"] = lon
            # record the velocity at the end of the current timestep (at the
            # point in time of the new trajectory point)
            ds_next_posn["u_traj"] = u_end
            ds_next_posn["v_traj"] = v_end
            points.append(ds_next_posn)

        if dir == "backward":
            # now we've integrated backwards we reverse the points
            points = points[::-1]

    ds_traj = xr.concat(points, dim="time").sortby("time")
    ds_traj["u_traj"].attrs = {"long_name": "zonal velocity", "units": "m/s"}
    ds_traj["v_traj"].attrs = {"long_name": "meridional velocity", "units": "m/s"}
    ds_traj["origin_lat"] = origin.lat
    ds_traj["origin_lon"] = origin.lon
    ds_traj["origin_datetime"] = origin.datetime
    return ds_traj


if __name__ == "__main__":
    main()
