import xarray as xr
import numpy as np


def extract_trajectory(lat0, lon0, t0, t_max, dt):
    times = [t0]
    while times[-1] < t_max:
        t_new = times[-1] + dt
        times.append(t_new)

    ds = xr.Dataset(coords=dict(time=times),)
    ds["lat"] = (("time",), lat0 * np.ones(len(times)))
    ds["lon"] = (("time",), lon0 * np.ones(len(times)))

    return ds
