from pathlib import Path

import xarray as xr
import numpy as np

from .geometry import longitude_set_meridian
from .constants import longitude_tolerance
from . import integrate_trajectory


def _era5_normalise_longitude(ds_to_normalise, ds_ref):
    """Normalise longitudes to be between 0 and 360 degrees
    This is needed because these are stored differently in the surface
    and model level data. Rounding up to 4 decimals seems to work for now,
    with more decimals misalignment has happenend. Would be good to sort
    out why this is the case.
    """
    ref_longitudes = longitude_set_meridian(ds_ref.coords["longitude"])
    these_longidues = longitude_set_meridian(ds_to_normalise.coords["longitude"])
    if len(ref_longitudes) == len(these_longidues):
        mean_square_error = np.mean(
            np.sqrt(
                (ref_longitudes - these_longidues) * (ref_longitudes - these_longidues)
            )
        )
        if mean_square_error < longitude_tolerance:
            out_longitudes = ref_longitudes
        else:
            out_longitudes = these_longidues
    else:
        out_longitudes = these_longidues
    ds_to_normalise.coords["longitude"] = (
        "longitude",
        out_longitudes,
        ds_to_normalise.coords["longitude"].attrs,
    )
    return ds_to_normalise


def main(data_path):
    """Dummy implementations for trajectory tool"""
    print(list(Path(data_path).glob("*")))

    files_model_an = f"{data_path}/an_model_*.nc"
    files_single_an = f"{data_path}/an_single_*.nc"
    files_model_fc = f"{data_path}/fc_model_*.nc"
    files_single_fc = f"{data_path}/fc_single_*.nc"

    ds_model_an = xr.open_mfdataset(files_model_an, combine="by_coords")
    # z needs to be dropped to prevent duplicity, lnsp is simply redundant
    ds_model_an = ds_model_an.drop_vars(["z", "lnsp"])
    ds_single_an = xr.open_mfdataset(files_single_an, combine="by_coords")
    ds_model_fc = xr.open_mfdataset(files_model_fc, combine="by_coords")
    ds_single_fc = xr.open_mfdataset(files_single_fc, combine="by_coords")
    ds_list = [ds_model_an, ds_single_an, ds_model_fc, ds_single_fc]
    for this_ds in ds_list:
        _era5_normalise_longitude(this_ds, ds_model_an)
    ds_merged = xr.merge(ds_list)

    dummy_trajectory_dict = {
        "lat_origin": 13.3,
        "lon_origin": -57.717,
        "datetime_origin": "2020-02-03T12:30",
        "backward_duration_hours": 3,
        "forward_duration_hours": 1,
        "nr_iterations_traj": 10,
        "velocity_strategy": "lower_troposphere_humidity_weighted",
        # "velocity_strategy": "prescribed_velocity",
        # "u_traj" : -6.0,
        # "v_traj" : -0.25,
        # "velocity_strategy": "velocity_at_height",
        # "velocity_height": 1000.0,
        "pres_cutoff_start": 60000.0,
        "pres_cutoff_end": 50000.0,
    }
    integrate_trajectory(ds_merged, dummy_trajectory_dict)

    # dummy_forcings_dict = {
        # "gradients_strategy": "both",
        # "mask": "ocean",
        # "traj_file": "ds_traj.nc",
        # "averaging_width": 4.0,
        # "w_cutoff_start": 70000.0,
        # "w_cutoff_end": 40000.0,
    # }
    # dummy_forcings(ds_merged, dummy_forcings_dict)


import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('data_path')
args = argparser.parse_args()

main(args.data_path)
