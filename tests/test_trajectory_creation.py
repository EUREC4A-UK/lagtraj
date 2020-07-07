from pathlib import Path
import dateutil.parser
import datetime
import numpy as np
import xarray as xr

import lagtraj.produce.lagrangian_trajectory
from lagtraj.utils import validation

import lagtraj.trajectory.create


def test_create_stationary_trajectory(ds_domain_test):
    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=40,
        lon=-10,
        datetime=datetime.datetime(year=2020, month=1, day=1, hour=10, minute=0),
    )

    da_times = ds_domain_test.time

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin, trajectory_type="eulerian", da_times=da_times,
    )

    validation.validate_trajectory(ds_traj)


def test_create_linear_trajectory(ds_domain_test):
    t0 = ds_domain_test.isel(time=0).time

    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=ds_domain_test.lat.mean(), lon=ds_domain_test.lon.mean(), datetime=t0,
    )

    da_times = ds_domain_test.time

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin, trajectory_type="linear", da_times=da_times, U=[5.0, 5.0,]
    )

    assert ds_traj.time.equals(da_times)

    validation.validate_trajectory(ds_traj)


def test_integrated_linear_trajectory(ds_domain_test):
    t0 = ds_domain_test.sel(time=slice("2020-01-01T00:00:00", None)).isel(time=0).time

    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=ds_domain_test.lat.mean(), lon=ds_domain_test.lon.mean(), datetime=t0,
    )

    t_min = t0.values
    dt = np.timedelta64("10", "m")
    t_max = t_min + dt * 3.0
    t_ = np.arange(t_min, t_max, dt)
    da_times = xr.DataArray(t_, dims=("time",), coords=dict(time=t_))

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin,
        trajectory_type="integrated",
        velocity_method="single_height_level",
        velocity_method_kwargs=dict(height=700.0, time_space_interpolation="nearest",),
        da_times=da_times,
        ds_domain=ds_domain_test,
    )

    assert ds_traj.time.equals(da_times)

    validation.validate_trajectory(ds_traj)


# testing command line interface
def test_cli():
    fn_out = "stationary_trajectory_test.nc"
    lagtraj.produce.lagrangian_trajectory.main(
        lat0=-10,
        lon0=40,
        t0=dateutil.parser.parse("2020-01-22T12:00"),
        t_max=dateutil.parser.parse("2020-01-24T12:00"),
        dt=datetime.timedelta(hours=4),
        trajectory_type="stationary",
        out_filename=fn_out,
    )

    p = Path(fn_out)
    assert p.exists()
    p.unlink()
