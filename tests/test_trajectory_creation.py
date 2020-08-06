import datetime

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
    ds_traj.attrs["name"] = "test_trajectory"
    ds_traj.attrs["domain_name"] = "test_domain_data"

    validation.validate_trajectory(ds_traj)


def test_create_linear_trajectory(ds_domain_test, ds_trajectory_linear):
    da_times = ds_domain_test.time
    ds_traj = ds_trajectory_linear
    assert ds_traj.time.equals(da_times)
    ds_traj.attrs["name"] = "test_trajectory"
    ds_traj.attrs["domain_name"] = "test_domain_data"
    validation.validate_trajectory(ds_traj)


def test_create_lagrangian_trajectory(ds_domain_test):
    da_times = ds_domain_test.time.isel(time=slice(-10, -5))
    t0 = da_times.isel(time=0)

    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=ds_domain_test.lat.mean(), lon=ds_domain_test.lon.mean(), datetime=t0,
    )

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin,
        trajectory_type="lagrangian",
        velocity_method="single_height_level",
        velocity_method_kwargs=dict(height=700.0,),
        da_times=da_times,
        ds_domain=ds_domain_test,
    )

    assert ds_traj.time.equals(da_times)
    ds_traj.attrs["name"] = "test_trajectory"
    ds_traj.attrs["domain_name"] = "test_domain_data"

    validation.validate_trajectory(ds_traj)
