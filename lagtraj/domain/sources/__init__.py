"""
Wrapper interface for dispatching calculations of auxiliary variables to the
correct functions depending on the data source a given dataset originated from.
"""
from .era5.interpolation import interpolate_to_height_levels as era5_hl_interp


def interpolate_to_height_levels(ds, height):
    """
    Some source data will not be define on "height levels" (i.e. the vertical
    coordinate represents values at the same height in meters), but instead
    might use a hybrid or pressure grid. This function calls the relevant
    interpolation routines to ensure the domain data exists on height levels.

    `height` is assumed to be in meters
    """
    data_source = ds.attrs.get("data_source")
    if data_source is None:
        raise Exception("Please define an attribute `data_source` on your domain data dataset"
                        " so that the correct method for transforming to height levels can"
                        " be used (e.g. `ds.attrs['data_source'] = 'era5')`")
    if data_source == "era5":
        return era5_hl_interp(ds_model_levels=ds, height=height)
    else:
        raise NotImplementedError(f"No method to inpolate domain data to height levels for"
                                  " `{data_source}` has been implemented")
