import datetime

import xarray as xr
import numpy as np

from ..utils import validation


def _extract_era5_forcing_profiles(ds_traj, required_variables, da_levels):
    if not da_levels.units == "Pa":
        raise NotImplementedError(
            "ERA5 data can currently only be generated on pressure levels"
        )

    ds_forcing_profiles = xr.Dataset(
        coords=dict(
            time=ds_traj.time,
            lat=ds_traj.lat,
            lon=ds_traj.lon,
            level=da_levels
        ),
        attrs=dict(description="ERA5 forcing data generated {}".format(
            datetime.datetime.now().isoformat()
        ))
    )

    ds_start = ds_traj.isel(time=0)
    t0, lat0, lon0 = ds_start.time, ds_start.lat, ds_start.lon

    # TODO: implement extraction from ERA5 data

    # extract trajectory data
    # for t in ds_traj.time ...

    raise NotImplementedError()


def _generate_dummy_forcing_profiles(ds_traj, required_variables, da_levels):
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


def main(fn_trajectory, source_data, out_filename):
    ds_traj = xr.open_dataset(fn_trajectory)
    validation.validate_trajectory(ds_traj)

    # TODO: move this to a yaml file
    da_levels = xr.DataArray(
        np.array([100e3, 90e3, 80e3, 70e3]),
        attrs=dict(long_name='pressure levels', units='Pa')
    )

    if source_data == 'ERA5':
        ds_forcing_profiles = _extract_era5_forcing_profiles(
            ds_traj=ds_traj, da_levels=da_levels
        )
    elif source_data == 'dummy':
        ds_forcing_profiles = _generate_dummy_forcing_profiles(
            ds_traj=ds_traj, required_variables=["ddt__qv"],
            da_levels=da_levels
        )
    else:
        raise NotImplementedError(source_data)

    validation.validate_forcing_profiles(ds_forcing_profiles)

    ds_forcing_profiles.to_netcdf(out_filename)
    print("Wrote forcing profiles to `{}`".format(out_filename))
