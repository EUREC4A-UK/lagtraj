import xarray as xr
from ..domain.sources.era5.load import ERA5DataSet


def append_dictionary_to_attrs(input_dictionary, output_ds, init_str=""):
    for k, v in input_dictionary.items():
        if isinstance(v, dict):
            append_dictionary_to_attrs(v, output_ds, init_str=str(k) + "_")
        elif isinstance(v, xr.Dataset) or isinstance(v, ERA5DataSet):
            output_ds.attrs[init_str + k] = str(type(v))
        elif isinstance(v, xr.DataArray):
            output_ds.attrs[init_str + k] = v.values
        else:
            output_ds.attrs[init_str + k] = v
