import datetime

import xarray as xr
import numpy as np

from ...utils import validation
from .dummy import extract_forcing_profiles as dummy_extract_forcing_profiles
from .era5 import extract_forcing_profiles as era5_extract_forcing_profiles


def main(fn_trajectory, source_data, out_filename):
    ds_traj = xr.open_dataset(fn_trajectory)
    validation.validate_trajectory(ds_traj)

    # TODO: move this to a yaml file
    da_levels = xr.DataArray(
        np.array([100e3, 90e3, 80e3, 70e3]),
        attrs=dict(long_name="pressure levels", units="Pa"),
    )

    if source_data == "ERA5":
        ds_forcing_profiles = era5_extract_forcing_profiles(
            ds_traj=ds_traj, required_variables=["ddt__qv"], da_levels=da_levels
        )
    elif source_data == "dummy":
        ds_forcing_profiles = dummy_extract_forcing_profiles(
            ds_traj=ds_traj, required_variables=["ddt__qv"], da_levels=da_levels
        )
    else:
        raise NotImplementedError(source_data)

    validation.validate_forcing_profiles(ds_forcing_profiles)

    ds_forcing_profiles.to_netcdf(out_filename)
    print("Wrote forcing profiles to `{}`".format(out_filename))
