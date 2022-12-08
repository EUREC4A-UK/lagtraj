import datetime
import functools
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from ....utils.units import round_time
from .. import MissingDomainData
from . import FILENAME_FORMAT, VERSION_FILENAME
from .utils import add_era5_global_attributes

MODEL_RUN_TYPES = ["an", "fc"]  # analysis and forecast runs
LEVEL_TYPES = ["model", "single"]  # need model and surface data


def _create_normalised_longitude(da_lon):
    """Normalise longitudes to be between 0 and 360 degrees
    This is needed because these are stored differently in the surface
    and model level data. Rounding up to 4 decimals seems to work for now,
    with more decimals misalignment has happenend. Would be good to sort
    out why this is the case.
    """

    def longitude_set_meridian(longitude):
        """Sets longitude to be between -180 and 180 degrees"""
        return (longitude + 180.0) % 360.0 - 180.0

    longitude_values = da_lon.data
    longitude_values_rounded_normed = np.round(
        longitude_set_meridian(longitude_values), decimals=4
    )
    return xr.DataArray(
        longitude_values_rounded_normed,
        dims=("longitude"),
        name="longitude",
        attrs=da_lon.attrs,
    )


def _load_single_file_preprocess(ds):
    """
    Preprocess a single file ensuring the time coordinate is monotonically
    increasing as expected. This sort is needed with some model-level files
    returned from ECMWF CDS as they have recently not been in correct time
    order (see discussion on https://github.com/EUREC4A-UK/lagtraj/issues/183
    for details)
    """
    time_index = ds._indexes.get("time").to_pandas_index()

    if not time_index.is_monotonic_increasing:
        ds = ds.sortby("time")
        warnings.warn(
            "loaded ERA5 file `{}` was sorted by the `time` coordinate during load "
            "to handle out-of-order data in this file"
        )

    return ds


def _find_datasets(data_path):
    datasets = {}

    for model_run_type in MODEL_RUN_TYPES:
        for level_type in LEVEL_TYPES:
            filename_format = FILENAME_FORMAT.format(
                model_run_type=model_run_type, level_type=level_type, date="*"
            )

            files = list(data_path.glob(filename_format))

            if len(files) == 0:
                raise MissingDomainData(
                    f"No files for {model_run_type} model run {level_type} "
                    f"level were found in {data_path}."
                )

            ds_ = xr.open_mfdataset(
                files, combine="by_coords", preprocess=_load_single_file_preprocess
            )
            # z needs to be dropped to prevent duplicity, lnsp is simply
            # redundant
            if model_run_type == "an" and level_type == "model":
                ds_ = ds_.drop_vars(["lnsp", "z"])
            da_lon = ds_.coords["longitude"]
            da_lon_normalised = _create_normalised_longitude(da_lon=da_lon)
            ds_ = ds_.assign_coords(dict(longitude=da_lon_normalised))
            ds_ = ds_.rename(dict(latitude="lat", longitude="lon"))

            dataset_identifier = f"{model_run_type}__{level_type}"
            datasets[dataset_identifier] = ds_
    return datasets


def _calc_creation_timestamp(data_path):
    """
    ERA5 datasets consist of multiple files and so working out a single
    creation time is not obvious. We need a datetime that is fixed once a
    dataset has been downloaded however. To ensure that the time doesn't change
    as the individual files are downloading we use the creation time of the
    file that was created first.
    """
    creation_times = []
    for model_run_type in MODEL_RUN_TYPES:
        for level_type in LEVEL_TYPES:
            filename_format = FILENAME_FORMAT.format(
                model_run_type=model_run_type, level_type=level_type, date="*"
            )
            files = data_path.glob(filename_format)

            creation_times += [fn.stat().st_ctime for fn in files]

    t = datetime.datetime.fromtimestamp(min(creation_times))

    # round to nearest second
    return round_time(t, num_seconds=1)


class ERA5DataSet(object):
    """
    Mimicks a xarray.Dataset containing files from ERA5, but without merging
    the individual files. Instead the files are sliced separately before
    merging each time a specific operation (for exameple interpolation) is
    applied
    """

    def __init__(self, data_path, selected_vars=[], datasets=[]):
        self.data_path = data_path
        self.datasets = datasets or _find_datasets(data_path=data_path)
        self._selected_vars = selected_vars
        self.attrs = {}

    def __getattr__(self, v):
        if v in self.attrs:
            return self.attrs[v]
        else:
            raise AttributeError(v)

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
            da_ = ds[v].rename("_temp")
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

    def __getitem__(self, item):
        requested_vars = type(item) == set and item or set(item)
        available_vars = self._selected_vars or self.data_vars
        missing_vars = requested_vars.difference(available_vars)
        if len(missing_vars) > 0:
            s = ", ".join(missing_vars)
            s2 = ", ".join(available_vars)
            raise Exception(
                f"Some of the variables you have requested ({s})"
                " aren't available in this dataset, the ones available"
                f" are: {s2}"
            )
        return ERA5DataSet(
            data_path=self.data_path,
            datasets=self.datasets,
            selected_vars=requested_vars,
        )

    @property
    @functools.lru_cache(maxsize=1)
    def time(self):
        return self._extra_var(v="time")

    @property
    @functools.lru_cache(maxsize=1)
    def lon(self):
        return self._extra_var(v="lon")

    @property
    @functools.lru_cache(maxsize=1)
    def lat(self):
        return self._extra_var(v="lat")

    @property
    @functools.lru_cache(maxsize=1)
    def data_vars(self):
        v = set()
        for ds in self.datasets.values():
            v = v.union(ds.data_vars)
        return v

    def sel(
        self, indexer=None, method=None, tolerance=None, drop=False, **indexers_kwargs
    ):
        if len(self._selected_vars) == 0:
            warnings.warn(
                "You are doing a selection on *all* variables available in this era5 dataset"
                " which is expensive (as they are loaded from many individual files). Instead"
                " select the variables you need (e.g. with `ds[('u', 'v')]`) before calling"
                " .sel to make a coordinated-based selection"
            )
        requested_variables = self._selected_vars or self.data_vars
        indexers_dims = indexers_kwargs.keys()
        datasets_slices = []

        for ds in self.datasets.values():
            variables = list(set(requested_variables).intersection(ds.data_vars))
            if len(variables) == 0:
                continue
            slices = {}
            dims = set(indexers_dims).intersection(ds.dims)
            for d in dims:
                s = indexers_kwargs[d]
                # ensure that slices always work if defined monotonically
                if ds[d][0] > ds[d][-1]:
                    slices[d] = slice(s.stop, s.start, s.step)
                else:
                    slices[d] = s

            ds_v_slice = ds[variables].sel(
                **slices,
                method=method,
                tolerance=tolerance,
                drop=drop,
            )
            ds_v_slice.load()
            # Change the long name of these variables, so the units are no
            # longer in contradicting with it (with cfchecker convections)
            for variable in ["sshf", "slhf", "sshf_local", "slhf_local"]:
                if variable in variables:
                    ds_v_slice[variable].attrs = {
                        "long_name": ds_v_slice[variable].long_name + " time integral",
                        "units": ds_v_slice[variable].units,
                    }
            datasets_slices.append(ds_v_slice)

        ds_ = xr.merge(datasets_slices, compat="override").load()
        ds_.attrs["data_source"] = "era5"
        return ds_

    def interp(self, kwargs, method="linear", **interp_to):
        """
        Implements xarray.interp by first slicing out the necessary data from
        individual era5 files, merging these and then interpolating
        """
        if len(self._selected_vars) == 0:
            warnings.warn(
                "You are doing an interpolation on *all* variables available in this era5 dataset"
                " which is expensive (as they are loaded from many individual files). Instead"
                " select the variables you need (e.g. with `ds[('u', 'v')]`) before calling"
                " .sel to make a coordinated-based selection"
            )
        requested_variables = self._selected_vars or self.data_vars

        if method not in ["linear", "nearest"]:
            raise NotImplementedError(
                "The interpolation method doesn't support"
                f" the {method} interpolation method,"
                " as not enough values are extracted during"
                " slicing."
            )

        interp_dims = interp_to.keys()
        datasets_slices = []
        for ds in self.datasets.values():
            # Note this needs to be a list, in order to refer to it when selecting from ds
            variables = list(set(requested_variables).intersection(ds.data_vars))
            if len(variables) == 0:
                continue
            slices = {}
            dims = set(interp_dims).intersection(ds.dims)
            for d in dims:
                # if a particular value is in a coordinate (for example
                # time) we just take that out, no need to make a slice
                d_vals_array = ds[d].values
                if np.array(interp_to[d]) in d_vals_array:
                    slices[d] = interp_to[d]
                    continue
                if type(interp_to[d]) == xr.core.dataarray.DataArray:
                    d_interp_val = interp_to[d].values
                else:
                    d_interp_val = interp_to[d]
                d_vals_smaller = d_vals_array[d_vals_array < d_interp_val]
                d_vals_greater = d_vals_array[d_vals_array > d_interp_val]
                at_left_edge = len(d_vals_smaller) == 0
                at_right_edge = len(d_vals_greater) == 0
                if at_left_edge or at_right_edge:
                    raise Exception(
                        "Requested interpolation at edge of domain"
                        f" (trying to access {d}={interp_to[d]}"
                        f" between {ds[d].min().values} and"
                        f" {ds[d].max().values}"
                    )
                d_slice_min = np.nanmax(d_vals_smaller)
                d_slice_max = np.nanmin(d_vals_greater)
                # Slice direction depends on how variables are ordered
                if d_vals_array[0] < d_vals_array[1]:
                    slices[d] = slice(d_slice_min, d_slice_max)
                else:
                    slices[d] = slice(d_slice_max, d_slice_min)
            ds_v_slice = ds[variables].sel(**slices)
            ds_v_slice.load()
            datasets_slices.append(ds_v_slice)

        ds_slice = xr.merge(datasets_slices, compat="override")

        # remove coords that we've picked already (matching exact values)
        extra_dims = list(set(interp_to.keys()).difference(ds_slice.dims))
        for d in extra_dims:
            del interp_to[d]
        ds_ = ds_slice.interp(**interp_to, kwargs=kwargs, method=method).load()
        ds_.attrs["data_source"] = "era5"
        return ds_


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
    creation_datetime = _calc_creation_timestamp(data_path)
    add_era5_global_attributes(ds, creation_datetime=creation_datetime)

    version_filename = Path(data_path) / VERSION_FILENAME
    if version_filename.exists():
        version = open(version_filename).read().strip()
    else:
        version = "unversioned"
    ds.attrs["version"] = version

    return ds
