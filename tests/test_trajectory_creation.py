from pathlib import Path
import dateutil.parser
import datetime

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
    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=40,
        lon=-10,
        datetime=datetime.datetime(year=2020, month=1, day=1, hour=10, minute=0),
    )

    da_times = ds_domain_test.time

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin, trajectory_type="linear", da_times=da_times, U=[5.0, 5.0,]
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
