"""
Tests for checking different options related to nudging when doing forcing
conversion to 

In this file we will be doing two things:

    1) Checking that the forcing step runs without error for different valid
    sets of parameters in the yaml-files that describe the conversion step
    2) When converting to DEPHY format we also check the output xr.Dataset
    (and so netCDF file), based on the conversion parameters used, whether the
    nudging related attributes have the correct values and the correct nudging
    variables are present (and have the correct values)
"""
import uuid
from pathlib import Path

import pytest
import xarray as xr
import yaml

import lagtraj
from lagtraj.forcings import ForcingLevelsDefinition, ForcingSamplingDefinition

FORCING_CONVERSION_BASE = """
levels_method   : copy
export_format   : kpt
comment         : Forcing and initial conditions for Lagrangian case
campaign        : EUREC4A
source_domain   : n/a
reference       : n/a
author          : s.j. boeing
modifications   : First version
case            : n/a
adv_temp        : 1
adv_theta       : 1
adv_thetal      : 1
adv_qv          : 1
adv_qt          : 1
adv_rv          : 1
adv_rt          : 1
rad_temp        : 0
rad_theta       : 0
rad_thetal      : 0
forc_omega      : 0
forc_w          : 1
forc_geo        : 1
surfaceType     : ocean
surfaceForcing  : ts
surfaceForcingWind               : z0_traj
"""


# 1. global attributes to check for
dict(
    # 0: off, -1: profile -2: run-time inversion height,
    # both constant and "fixed_height" us the "profile" option (-1), and we
    # will introduce our own option (-2) for "runtime_inversion_height"
    nudging_u,
    nudging_v,
    nudging_temp,
    nudging_theta,
    nudging_thetal,
    nudging_qv,
    nudging_qt,
    nudging_rv,
    nudging_rt,
)

# 2. variables to look for within the dataset
# shape: [time, height, lat, lon]

variables = """
nudging_inv_u_traj
nudging_inv_v_traj
nudging_inv_temp_traj
nudging_inv_theta_traj
nudging_inv_thetal_traj
nudging_inv_qv_traj
nudging_inv_qt_traj
nudging_inv_rv_traj
nudging_inv_rt_traj
"""


nudging_method_scalars: constant
nudging_method_scalars: runtime_inversion_height
nudging_method_scalars: fixed_height

if method == "constant":
    raise ValueError("not implemented")
    # all values for momentum/scalars should be the inverse of the timescale
    # given by nudging_timescale_scalars/nudging_timescale_momentum
elif method == "fixed_height":
    raise ValueError("not implemented")
    # below the fixed height the nudging_inv_* values are zero and at the very
    # domain top as method=="constant" above
elif method == "runtime_inversion_height":
    raise ValueError("not implemented")
    # check that variables for nudging profiles aren't present
elif method == "off":
    raise ValueError("not implemented")
    # no variable presents


# 3. ensure parameters in yaml file describing conversion are present as
# attribute in output netCDF file


VALID_YAML_EXAMPLES = [
    """
nudging_method_scalars: off
""",
    """
nudging_method_scalars: constant
nudging_timescale_scalars: 10
""",
    """
nudging_method_scalars: runtime_inversion_height
nudging_timescale_scalars: 10800
nudging_transition_thickness_scalars: 500.0
nudging_shape_scalars: cos
""",
    """
nudging_method_scalars: fixed_height
nudging_timescale_scalars: 10800
nudging_transitions_thickness_scalars: 500.0
nudging_shape_scalars: cos
nudging_above_height_scalars: 1000.0
""",
]


@pytest.fixture(scope="session")
def example_forcing(ds_domain_test, ds_trajectory_linear):
    fp = Path("/tmp/lolforcing.nc")
    if fp.exists():
        return xr.open_dataset(fp)

    # just use five timesteps to make the tests execute faster
    ds_traj = ds_trajectory_linear.isel(time=slice(0, 2))
    ds_domain = ds_domain_test

    levels_definition = ForcingLevelsDefinition(
        method="exponential",
        n_levels=10,
        z_top=6.0e3,
        dz_min=100.0,
    )

    sampling_definition = ForcingSamplingDefinition(
        gradient_method="boundary",
        advection_velocity_sampling_method="domain_mean",
        averaging_width=2.0,
        time_sampling_method="domain_data",
        mask="ocean_only",
    )

    forcing_name = "test_forcing"

    ds_forcing = lagtraj.forcings.create.make_forcing(
        ds_trajectory=ds_traj,
        ds_domain=ds_domain,
        levels_definition=levels_definition,
        sampling_method=sampling_definition,
    )

    # to enable export to file of a forcing we much give it a name
    ds_forcing.attrs["name"] = forcing_name

    ds_forcing.to_netcdf(fp)

    return ds_forcing


def _store_forcing_conversion_input_definition(
    conversion_yaml, root_data_path, forcing_name
):
    """
    Store the yaml content that describes the forcing conversion into the
    location where it can be used to apply a forcing conversion
    """
    # use a random unique identifier to name this conversion
    example_identifier = uuid.uuid4().hex
    conversion_name = f"forcing_nudging__{example_identifier}"

    # save the forcing conversion yaml-file to the test-data directory so that
    # we can use it below to convert the forcing
    conversion_params = yaml.load(conversion_yaml, Loader=yaml.FullLoader)
    conversion_filename = f"{forcing_name}.{conversion_name}.yaml"
    conversion_fp = Path(root_data_path) / "forcings" / conversion_filename
    conversion_fp.parent.mkdir(exist_ok=True, parents=True)
    with open(conversion_fp, "w", encoding="utf-8") as fh:
        yaml.dump(conversion_params, fh)

    return conversion_name


@pytest.mark.parametrize("example_yaml", VALID_YAML_EXAMPLES)  # noqa
def test_forcing_conversion_config(example_forcing, example_yaml):
    # we export to a temporary directory so that we don't clobber local
    # files or have parallel tests clobber each other
    ds_forcing = example_forcing

    root_data_path = Path("/tmp/leiftest")

    # with tempfile.TemporaryDirectory() as root_data_path:
    conversion_yaml = FORCING_CONVERSION_BASE + example_yaml
    conversion_name = _store_forcing_conversion_input_definition(
        conversion_yaml=conversion_yaml,
        forcing_name=ds_forcing.name,
        root_data_path=root_data_path,
    )

    def run_conversion():
        lagtraj.forcings.conversion.process.export_for_target(
            ds_forcing=ds_forcing,
            conversion_name=conversion_name,
            root_data_path=root_data_path,
        )

    run_conversion()
