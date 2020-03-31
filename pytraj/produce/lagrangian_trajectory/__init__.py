import datetime

import xarray as xr
import numpy as np
import dateutil.parser

from ...utils import validation
from .era5 import extract_trajectory as extract_era5_trajectory
from .stationary import extract_trajectory as extract_stationary_trajectory


def main(lat0, lon0, t0, t_max, dt, trajectory_type, out_filename):
    if trajectory_type == "stationary":
        ds = extract_stationary_trajectory(
            lat0=lat0, lon0=lon0, t_max=t_max, t0=t0, dt=dt
        )
    elif trajectory_type == "ERA5":
        ds = extract_era5_trajectory(
            lat0=lat0, lon0=lon0, t_max=t_max, t0=t0, dt=dt
        )
    else:
        raise NotImplementedError(trajectory_type)

    validation.validate_trajectory(ds)
    ds.to_netcdf(out_filename)
    print("Wrote trajectory to {}".format(out_filename))
