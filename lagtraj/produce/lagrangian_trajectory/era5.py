import datetime

import xarray as xr
import numpy as np
import dateutil.parser


def extract_trajectory(lat0, lon0, t0, t_max, dt):
    if dt != datetime.timedelta(hours=3):
        raise NotImplementedError("ERA5 trajectories currently have to have"
                                  " 3-hour timestep")

    # implement fetching data and integrating trajectory
    raise NotImplementedError("ERA5 trajectory integration not yet completed")
