import xarray as xr
import numpy as np

from . import FILENAME_FORMAT


def _era_5_normalise_longitude(ds):
    """Normalise longitudes to be between 0 and 360 degrees
    This is needed because these are stored differently in the surface
    and model level data. Rounding up to 4 decimals seems to work for now,
    with more decimals misalignment has happenend. Would be good to sort
    out why this is the case.
    """

    def longitude_set_meridian(longitude):
        """Sets longitude to be between -180 and 180 degrees"""
        return (longitude + 180.0) % 360.0 - 180.0

    ds.coords["longitude"] = (
        "longitude",
        np.round(longitude_set_meridian(ds.coords["longitude"]), decimals=4),
        ds.coords["longitude"].attrs,
    )
    return ds


def load_data(data_path):
    datasets = {}

    model_run_types = ["an", "fc"]  # analysis and forecast runs
    level_types = ["model", "single"]  # need model and surface data

    for model_run_type in model_run_types:
        datasets_run = []
        for level_type in level_types:
            filename_format = FILENAME_FORMAT.format(
                model_run_type=model_run_type, level_type=level_type, date="*"
            )

            files = data_path.glob(filename_format)

            ds_ = xr.open_mfdataset(files, combine="by_coords")
            # z needs to be dropped to prevent duplicity, lnsp is simply
            # redundant
            if model_run_type == "an" and level_type == "model":
                ds_ = ds_.drop_vars(["lnsp"])
            ds_ = _era_5_normalise_longitude(ds=ds_)
            datasets_run.append(ds_)
        ds_run = xr.merge(datasets_run, compat="override")
        datasets[model_run_type] = ds_run

    ds = xr.merge(datasets.values(), compat="override")
    ds = ds.rename(dict(latitude="lat", longitude="lon"))

    return ds
