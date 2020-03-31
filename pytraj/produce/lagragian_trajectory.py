import datetime

import xarray as xr
import numpy as np
import dateutil.parser


def _make_stationary_trajectory(lat0, lon0, t0, t_max, dt):
    times = [t0]
    while times[-1] < t_max:
        t_new = times[-1] + dt
        times.append(t_new)

    ds = xr.Dataset(
        coords=dict(time=times),
    )
    ds['lat'] = (('time',), lat0*np.ones(len(times)))
    ds['lon'] = (('time',), lon0*np.ones(len(times)))

    return ds


def _extract_trajectory_from_era5(lat0, lon0, t0, t_max, dt):
    # TODO: move this to a separate source file once a lot of functionality is
    # added
    if dt != datetime.timedelta(hours=3):
        raise NotImplementedError("ERA5 trajectories currently have to have"
                                  " 3-hour timestep")

    # implement fetching data and integrating trajectory
    raise NotImplementedError("ERA5 trajectory integration not yet completed")


def main(lat0, lon0, t0, t_max, dt, trajectory_type, out_filename):
    if trajectory_type == "stationary":
        ds = _make_stationary_trajectory(
            lat0=lat0, lon0=lon0, t_max=t_max, t0=t0, dt=dt
        )
    elif trajectory_type == "ERA5":
        ds = _extract_trajectory_from_era5(
            lat0=lat0, lon0=lon0, t_max=t_max, t0=t0, dt=dt
        )
    else:
        raise NotImplementedError(trajectory_type)

    ds.to_netcdf(out_filename)
    print("Wrote trajectory to {}".format(out_filename))


if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('lat0', help='starting latitude')
    argparser.add_argument('lon0', help='starting longitude')
    argparser.add_argument('t0', help='starting time',
                           type=dateutil.parser.parse)
    argparser.add_argument('t_max', help='end time',
                           type=dateutil.parser.parse)
    argparser.add_argument('dt', help='time increment')
    argparser.add_argument('--out-filename', help='output filename (netCDF)')

    args = argparser.parse_args()

    main(**dict(args))

