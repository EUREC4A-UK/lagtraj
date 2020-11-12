import pytest


import lagtraj.forcings.create
from lagtraj.forcings import ForcingLevelsDefinition, ForcingSamplingDefinition
from lagtraj.utils import validation


@pytest.mark.parametrize(
    "gradient_method, advection_velocity_method",
    [
        ("boundary", "mean"),
        ("regression", "mean"),
        ("boundary", "local"),
        ("regression", "local"),
    ],
)
def test_create_forcing_linear_trajectory(
    ds_domain_test, ds_trajectory_linear, gradient_method, advection_velocity_method
):
    # just use five timesteps to make the tests execute faster
    ds_traj = ds_trajectory_linear.isel(time=slice(0, 5))
    ds_domain = ds_domain_test

    levels_definition = ForcingLevelsDefinition(
        method="exponential", n_levels=10, z_top=6.0e3, dz_min=100.0,
    )
    sampling_definition = ForcingSamplingDefinition(
        gradient_method=gradient_method,
        advection_velocity_method=advection_velocity_method,
        averaging_width=2.0,
        time_sampling_method="domain_data",
        mask="ocean_only",
    )
    ds_forcing = lagtraj.forcings.create.make_forcing(
        ds_trajectory=ds_traj,
        ds_domain=ds_domain,
        levels_definition=levels_definition,
        sampling_method=sampling_definition,
    )

    validation.validate_forcing_profiles(ds_forcing)
    validation.check_for_ncview_warnings(ds=ds_forcing)
