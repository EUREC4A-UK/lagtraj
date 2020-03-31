import datetime

import xarray as xr
import numpy as np


def extract_forcing_profiles(ds_traj, required_variables, da_levels):
    ds_start = ds_traj.isel(time=0)
    t0, lat0, lon0 = ds_start.time, ds_start.lat, ds_start.lon

    ds_forcing_profiles = xr.Dataset(
        coords=dict(
            time=ds_traj.time,
            lat=ds_traj.lat,
            lon=ds_traj.lon,
            level=da_levels
        ),
        attrs=dict(description="Dummy forcing data generated {}".format(
            datetime.datetime.now().isoformat()
        ))
    )

    for req_var in required_variables:
        vals = np.random.random(
            (len(ds_forcing_profiles.time), 
             len(ds_forcing_profiles.lat),
             len(ds_forcing_profiles.lon),
             len(ds_forcing_profiles.level))
        )

        ds_forcing_profiles[req_var] = xr.DataArray(
            data=vals,
            dims=('time', 'lat', 'lon', 'level'),
            attrs=dict(
                units="1",
                long_name='foobar values for {}'.format(req_var)
            )
        )

    return ds_forcing_profiles
