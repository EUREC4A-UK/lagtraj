import datetime
import warnings
from pathlib import Path

import numpy as np
import tqdm
import xarray as xr

from .. import DEFAULT_ROOT_DATA_PATH
from ..domain.download import download_complete
from ..domain.load import load_data as load_domain_data
from ..utils import optional_debugging, validation
from ..utils.units import fix_units
from ..utils.xarray import create_attributes_dictionary
from . import build_data_path, extrapolation
from .load import load_definition

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
        ds_traj = create_eulerian_trajectory(origin=origin, da_times=da_times)
    elif trajectory_type == "linear":
        if "U" not in kwargs:
            raise Exception(
                "To use the `linear` trajectory integration you"
                " must provide a velocity `U`"
            )
        ds_traj = create_linear_trajectory(origin=origin, da_times=da_times, **kwargs)
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
        ds_traj = create_lagrangian_trajectory(
            origin=origin, da_times=da_times, **kwargs
        )
    else:
        raise NotImplementedError(f"`{trajectory_type}` trajectory type not available")

    # before we make the attributes rename domain so that it isn't called
    # `ds_domain` in the attributes
    if "ds_domain" in kwargs:
        kwargs["domain"] = kwargs.pop("ds_domain")
    ds_traj.attrs.update(
        create_attributes_dictionary(trajectory_type=trajectory_type, **kwargs)
    )
    return ds_traj


def main(data_path, trajectory_name):
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
    elif (
        traj_definition.type == "eulerian" or traj_definition.type == "linear"
    ) and traj_definition.timestep == "domain_data":
        ds_domain = load_domain_data(
            root_data_path=data_path, name=traj_definition.domain
        )
        kwargs["ds_domain"] = ds_domain
    else:
        ds_domain = None

    ds_trajectory = create_trajectory(**kwargs)

    # store all the parameters of how this trajectory was created into the
    # attributes of the dataset
    configuration = dict(name=trajectory_name)
    if ds_domain is not None:
        configuration["domain"] = ds_domain.attrs
    ds_trajectory.attrs.update(
        create_attributes_dictionary(traj_definition, **configuration)
    )

    trajectory_data_path = build_data_path(
        root_data_path=data_path, trajectory_name=traj_definition.name
    )
    validation.validate_trajectory(ds_traj=ds_trajectory)
    encoding = validation.build_valid_encoding(ds=ds_trajectory)
    ds_trajectory.to_netcdf(trajectory_data_path, encoding=encoding)
    print("Saved trajectory to `{}`".format(trajectory_data_path))


def _get_times_from_domain(trajectory_definition, root_data_path):
    """
    Return `xr.DataArray` of `datetime`s as a subset of the times at which the
    requested domain data is available, ensuring that the provided trajectory
    definition can be interpolated through its full duration from start to
    finish

    We need to be careful here to ensure here that the earliest and latest
    time that the trajectory was requested at is within the sampling times we
    actually use. For example, when using ERA5 domain data (which is available
    on the hour) while at the same setting a origin time for the trajectory
    15min past the hour and integrating +-2 hrs, we need to ensure that the
    domain data includes the time-steps for the hour immediately before and after
    the trajectory's full span in time, and so we need to find the latest time
    which is before the start of the trajectory and the earliest time that is
    after the trajectory end.

    hours past midnight:   0...1...2...3...4...5...6...7
    domain times:          x   x   x   x   x   x   x   x
    trajectory times:           s       o       e
    correct time-span:         |-------------------|
    """
    t0 = trajectory_definition.origin.datetime
    t_min_traj = t0 - trajectory_definition.duration.backward
    t_max_traj = t0 + trajectory_definition.duration.forward

    if not download_complete(
        root_data_path=root_data_path,
        domain_name=trajectory_definition.domain,
        start_date=t_min_traj.date(),
        end_date=t_max_traj.date(),
    ):
        warnings.warn(
            "Some of the data for the selected domain"
            f" ({trajectory_definition.domain}) hasn't been"
            " downloaded yet"
        )
    ds_domain = load_domain_data(
        root_data_path=root_data_path, name=trajectory_definition.domain
    )

    # have to convert to `np.datetime64` to be able to make comparisons below,
    # can't compare `datetime.datetime` with `xr.DataArray` of
    # `np.datetime64[ns]` objects
    if np.datetime64(t_min_traj) < ds_domain.time.min():
        raise Exception(
            "You selected to use the domain data for timesteps"
            " in the trajectory, but the start-time for the trajectory"
            " is beyond the beginning of which domain data is currently"
            " available for the `{trajectory_definition.domain}` domain"
        )

    if ds_domain.time.max() < np.datetime64(t_max_traj):
        raise Exception(
            "You selected to use the domain data for timesteps"
            " in the trajectory, but the end-time for the trajectory"
            " is beyond the end of which domain data is currently"
            " available for the `{trajectory_definition.domain}` domain"
        )

    da_tmin_domain = ds_domain.time.sel(time=slice(None, t_min_traj)).max()
    da_tmax_domain = ds_domain.time.sel(time=slice(t_max_traj, None)).min()
    da_times = ds_domain.time.sel(time=slice(da_tmin_domain, da_tmax_domain))
    if da_times.count() == 0:
        raise Exception(
            "You selected to use the domain data for timesteps"
            " in the trajectory, but in the time interval selected"
            " for the trajectory ({}, {}) there are is no domain data"
            " (time range: {} to {})".format(
                t_min_traj,
                t_max_traj,
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
    ds["origin_datetime"] = origin.datetime
    ds["origin_lat"] = xr.DataArray(lat0)
    ds["origin_lon"] = xr.DataArray(lon0)
    ds["lat"] = ("time"), ds.origin_lat.item() * np.ones(len(ds.time))
    ds["lon"] = ("time"), ds.origin_lon.item() * np.ones(len(ds.time))
    ds["u_traj"] = (
        ("time",),
        np.zeros(len(ds.time)),
        {"long_name": "zonal velocity", "units": "m s**-1"},
    )
    ds["v_traj"] = (
        ("time",),
        np.zeros(len(ds.time)),
        {"long_name": "meridional velocity", "units": "m s**-1"},
    )
    ds["origin_lat"].attrs = {
        "long_name": "latitude of trajectory origin",
        "units": "degrees_north",
        "info": "the reference point is the space-time coordinate from which the trajectory is calculated",
    }
    ds["origin_lon"].attrs = {
        "long_name": "longitude of trajectory origin",
        "units": "degrees_east",
        "info": "the reference point is the space-time coordinate from which the trajectory is calculated",
    }
    ds["origin_datetime"].attrs["long_name"] = "time of trajectory origin"
    ds["origin_datetime"].attrs[
        "info"
    ] = "the reference point is the space-time coordinate from which the trajectory is calculated"
    fix_units(ds)
    return ds


def create_linear_trajectory(origin, da_times, U):
    """Create linear trajectory from origin point using constant velocity"""

    def extrapolation_func(lat, lon, t0, dt):
        if dt > 0:
            s = 1.0
        else:
            s = -1.0

        lat_new, lon_new = extrapolation.extrapolate_posn_with_fixed_velocity(
            lat=lat,
            lon=lon,
            u_vel=s * U[0],
            v_vel=s * U[1],
            dt=s * dt,
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

        for t in tqdm.tqdm(da_integrate_times, desc=dir):
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
                points[0]["v_traj"] = v_start

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
    ds_traj["u_traj"].attrs = {"long_name": "zonal velocity", "units": "m s**-1"}
    ds_traj["v_traj"].attrs = {"long_name": "meridional velocity", "units": "m s**-1"}
    ds_traj["origin_lat"] = origin.lat
    ds_traj["origin_lon"] = origin.lon
    ds_traj["origin_datetime"] = origin.datetime
    ds_traj["lat"].attrs = {"long_name": "latitude", "units": "degrees_north"}
    ds_traj["lon"].attrs = {"long_name": "longitude", "units": "degrees_east"}
    ds_traj["origin_lat"].attrs = {
        "long_name": "latitude of trajectory origin",
        "units": "degrees_north",
        "info": "the reference point is the space-time coordinate from which the trajectory is calculated",
    }
    ds_traj["origin_lon"].attrs = {
        "long_name": "longitude of trajectory origin",
        "units": "degrees_east",
        "info": "the reference point is the space-time coordinate from which the trajectory is calculated",
    }
    ds_traj["origin_datetime"].attrs["long_name"] = "time of trajectory origin"
    ds_traj["origin_datetime"].attrs[
        "info"
    ] = "the reference point is the space-time coordinate from which the trajectory is calculated"
    fix_units(ds_traj)
    return ds_traj


def _make_cli_argparser():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("trajectory")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument("--debug", default=False, action="store_true")
    return argparser


def cli(args=None):
    """
    Function called with arguments passed from the command line when making
    trajectories through the CLI. When `args==None` they will be taken from
    `sys.argv`
    """
    argparser = _make_cli_argparser()
    args = argparser.parse_args(args=args)

    with optional_debugging(args.debug):
        main(data_path=args.data_path, trajectory_name=args.trajectory)


def has_data_for_cli_command(args):
    argparser = _make_cli_argparser()
    args = argparser.parse_args(args=args)

    data_path = args.data_path
    trajectory_name = args.trajectory
    return required_data_available(data_path=data_path, trajectory_name=trajectory_name)


def required_data_available(data_path, trajectory_name):
    traj_definition = load_definition(root_data_path=data_path, name=trajectory_name)

    if traj_definition.domain is not None:
        ds_domain = load_domain_data(
            root_data_path=data_path, name=traj_definition.domain
        )
        t0 = traj_definition.origin.datetime
        traj_t_min = t0 - traj_definition.duration.backward
        traj_t_max = t0 + traj_definition.duration.forward

        def dt64_to_dt(v):
            # to make a datetime.datetime we first have to remove the
            # nanosecond precision
            return v.astype("datetime64[s]").astype(datetime.datetime)

        data_t_min = dt64_to_dt(ds_domain.time.min().values)
        data_t_max = dt64_to_dt(ds_domain.time.max().values)

        has_time_range = data_t_min <= traj_t_min and traj_t_max <= data_t_max
        return has_time_range

    return True


if __name__ == "__main__":
    cli()
