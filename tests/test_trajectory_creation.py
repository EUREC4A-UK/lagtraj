import datetime
import inspect
from pathlib import Path

import isodate
import numpy as np
import xarray as xr
import yaml

import lagtraj.trajectory.create
from lagtraj.utils import validation


def test_create_stationary_trajectory(ds_domain_test):
    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=40,
        lon=-10,
        datetime=datetime.datetime(year=2020, month=1, day=1, hour=10, minute=0),
    )

    da_times = ds_domain_test.time

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin,
        trajectory_type="eulerian",
        da_times=da_times,
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
    validation.check_for_ncview_warnings(ds=ds_traj)


def test_create_lagrangian_trajectory(ds_domain_test):
    da_times = ds_domain_test.time.isel(time=slice(-10, -5))
    t0 = da_times.isel(time=0)

    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=ds_domain_test.lat.mean(),
        lon=ds_domain_test.lon.mean(),
        datetime=t0,
    )

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin,
        trajectory_type="lagrangian",
        velocity_method="single_height_level",
        velocity_method_kwargs=dict(
            height=700.0,
        ),
        da_times=da_times,
        ds_domain=ds_domain_test,
    )

    assert ds_traj.time.equals(da_times)
    ds_traj.attrs["name"] = "test_trajectory"
    ds_traj.attrs["domain_name"] = "test_domain_data"

    validation.validate_trajectory(ds_traj)


TIMESTEP_TEST_YAML_TEMPLATE = """
lat_origin: 13.0
lon_origin: -54.0
trajectory_type: eulerian
datetime_origin: 2020-02-02T12:15
forward_duration: PT3H
backward_duration: null
domain: lagtraj://{domain_name}
timestep: domain_data
"""


def test_check_for_all_timesteps(ds_domain_test):
    # name the trajectory yaml file using the current function's name
    fn_name = inspect.getframeinfo(inspect.currentframe()).function
    mod_name = __name__
    name = f"{mod_name}__{fn_name}"

    domain_name = ds_domain_test.name

    traj_example_yaml = TIMESTEP_TEST_YAML_TEMPLATE.format(domain_name=domain_name)
    traj_params = yaml.load(traj_example_yaml)

    data_path_root = Path(ds_domain_test.data_path).parent.parent

    traj_fp = Path(data_path_root) / "trajectories" / f"{name}.yaml"
    traj_fp.parent.mkdir(exist_ok=True, parents=True)
    with open(traj_fp, "w") as fh:
        yaml.dump(traj_params, fh)

    lagtraj.trajectory.create.main(data_path=data_path_root, trajectory_name=name)

    traj_fp_nc = traj_fp.parent / traj_fp.name.replace(".yaml", ".nc")
    ds_traj = xr.open_dataset(traj_fp_nc)

    t_origin = isodate.parse_datetime(traj_params["datetime_origin"])
    dt = isodate.parse_duration(traj_params["forward_duration"])
    t_max = t_origin + dt

    da_domain_times = ds_domain_test.time

    # we're using the timesteps of the source domain data, but the times that
    # result should contain the time span we have requested
    da_t_min_expected = ds_domain_test.sel(time=slice(None, t_origin)).time.max()
    da_t_max_expected = ds_domain_test.sel(time=slice(t_max, None)).time.min()

    da_t_expected = da_domain_times.sel(
        time=slice(da_t_min_expected, da_t_max_expected)
    )

    assert np.testing.assert_allclose(da_t_expected, ds_traj.time)
