import datetime

import xarray as xr
import numpy as np


def extract_forcing_profiles(ds_traj, required_variables, da_levels):
    if not da_levels.units == "Pa":
        raise NotImplementedError(
            "ERA5 data can currently only be generated on pressure levels"
        )

    ds_forcing_profiles = xr.Dataset(
        coords=dict(
            time=ds_traj.time, lat=ds_traj.lat, lon=ds_traj.lon, level=da_levels
        ),
        attrs=dict(
            description="ERA5 forcing data generated {}".format(
                datetime.datetime.now().isoformat()
            )
        ),
    )

    ds_start = ds_traj.isel(time=0)
    t0, lat0, lon0 = ds_start.time, ds_start.lat, ds_start.lon

    # TODO: implement extraction from ERA5 data

    # extract trajectory data
    # for t in ds_traj.time ...

    raise NotImplementedError()
