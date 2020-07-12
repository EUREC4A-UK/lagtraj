import xarray as xr
import numpy as np

import functools
import operator
import warnings

from . import FILENAME_FORMAT


MODEL_RUN_TYPES = ["an", "fc"]  # analysis and forecast runs
LEVEL_TYPES = ["model", "single"]  # need model and surface data


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


def _find_datasets(data_path):
    datasets = {}

    for model_run_type in MODEL_RUN_TYPES:
        for level_type in LEVEL_TYPES:
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
            ds_ = ds_.rename(dict(latitude="lat", longitude="lon"))

            dataset_identifier = f"{model_run_type}__{level_type}"
            datasets[dataset_identifier] = ds_
    return datasets


class ERA5DataSet(object):
    """
    Mimicks a xarray.Dataset containing files from ERA5, but without merging
    the individual files. Instead the files are sliced separately before
    merging each time a specific operation (for exameple interpolation) is
    applied
    """

    def __init__(self, data_path):
        self.data_path = data_path
        self.datasets = _find_datasets(data_path=data_path)

    def _extra_var(self, v):
        """
        Extract a particular variable across all files and extract the part
        where there is overlap between analysis and forecast files, used to
        extract a common time, lat and lon reference for all files
        """
        das_by_run_type = {}
        for ds_identifier, ds in self.datasets.items():
            model_run_type, _ = ds_identifier.split("__")
            # rename is required otherwise xarray complains when merging
            # because the dataarrays are promoted to datasets before merging,
            # and the dataset name would otherwise be the same as a dimension
            # (which is not allowed)
            da_ = ds[v].rename('_temp')
            das_by_run_type.setdefault(model_run_type, []).append(da_)

        # outer join first because era5 files are separated in time
        dss = []
        for model_run_type, das_run_type in das_by_run_type.items():
            ds_ = xr.merge(das_run_type, join="outer")
            dss.append(ds_)

        # finally, use inner-join to only return range for which all source
        # files have data
        da_combined = xr.merge(dss, join="inner")["_temp"]
        return da_combined

    @property
    def time(self):
        return self._extra_var(v="time")

    @property
    def lon(self):
        return self._extra_var(v="lon")

    @property
    def lat(self):
        return self._extra_var(v="lat")

    @property
    def data_vars(self):
        v = set()
        for ds in self.datasets.values():
            v = v.union(ds.data_vars)
        return v

    @property
    def sel(
        self, indexer=None, method=None, tolerance=None, drop=False, **indexers_kwargs
    ):
        warnings.warn(
            "Subsetting an entire era5 dataset is computationally expensive."
            " Consider accessing the the time, lat, lon attributes directly"
            " if these qre needed."
        )
        indexers_dims = indexers_kwargs.keys()
        datasets_slices = []
        for ds in self.datasets.values():
            das = []
            for v in ds.data_vars:
                da_v = ds[v]
                dims = set(indexers_dims).intersection(da_v.dims)

                slices = {}
                for d in dims:
                    slices[d] = indexers_kwargs[d]

                da_v_slice = da_v.sel(
                    **slices, method=method, tolerance=tolerance, drop=drop,
                )
                das.append(da_v_slice)

            ds_slice = xr.merge(das)
            datasets_slices.append(ds_slice)

        return xr.merge(datasets_slices, compat="override").compute()

    def interp(self, kwargs, **interp_to):
        """
        Implements xarray.interp by first slicing out the necessary data from
        individual era5 files, merging these and then interpolating
        """
        idx_padding = 1
        interp_dims = interp_to.keys()
        datasets_slices = []
        for ds in self.datasets.values():
            das = []
            for v in ds.data_vars:
                da_v = ds[v]
                dims = set(interp_dims).intersection(da_v.dims)

                slices = {}
                for d in dims:
                    # need to handle case where order isn't monototnically
                    # increasing
                    if da_v[d].values[0] > da_v[d].values[-1]:
                        dir = -1
                    else:
                        dir = 1
                    da_coord_left = da_v[d].sel(**{d: slice(None, interp_to[d], dir)})
                    da_coord_right = da_v[d].sel(**{d: slice(interp_to[d], None, dir)})

                    at_left_edge = da_coord_left.count() <= idx_padding
                    at_right_edge = da_coord_right.count() <= idx_padding
                    if at_left_edge or at_right_edge:
                        raise Exception(
                            "Requested interpolation at edge of domain"
                            f" (trying to access {d}={interp_to[d].values}"
                            f" between {da_v[d].min().values} and"
                            f" {da_v[d].max().values}"
                        )

                    d_edge_min = da_coord_left.isel(**{d: -1 - idx_padding})
                    d_edge_max = da_coord_right.isel(**{d: idx_padding})
                    if dir == -1:
                        slices[d] = slice(d_edge_max, d_edge_min)
                    else:
                        slices[d] = slice(d_edge_min, d_edge_max)
                    assert d_edge_min <= interp_to[d] <= d_edge_max
                    # make sure we've actually selected some data
                    assert functools.reduce(operator.mul, da_v.shape) > 0

                da_v_slice = da_v.sel(**slices)
                das.append(da_v_slice)

            ds_slice = xr.merge(das)
            datasets_slices.append(ds_slice)

        ds_slice = xr.merge(datasets_slices, compat="override").compute()
        return ds_slice.interp(**interp_to, kwargs=kwargs)


def _load_naive(data_path):
    """
    Load and merge all era5 files at once. NOTE: this uses a lot more memory
    and is slower than using the ERA5Dataset method
    """
    datasets = _find_datasets(data_path=data_path)

    ds = xr.merge(datasets.values(), compat="override")
    return ds


def load_data(data_path, use_lazy_loading=False):
    if use_lazy_loading:
        ds = _load_naive(data_path)
    else:
        ds = ERA5DataSet(data_path)

    return ds
