import xarray as xr
from ..domain.sources.era5.load import ERA5DataSet


def create_attributes_dictionary(input_dictionary, init_str=""):
    attributes_dictionary = {}
    for k, v in input_dictionary.items():
        if isinstance(v, dict):
            attributes_dictionary.update(
                create_attributes_dictionary(v, init_str=str(k) + "_")
            )
        elif isinstance(v, xr.Dataset) or isinstance(v, ERA5DataSet):
            attributes_dictionary[f"{init_str}{k}"] = v.__class__.__name__
        elif isinstance(v, xr.DataArray):
            attributes_dictionary[f"{init_str}{k}"] = v.values
        else:
            attributes_dictionary[f"{init_str}{k}"] = v
    return attributes_dictionary
