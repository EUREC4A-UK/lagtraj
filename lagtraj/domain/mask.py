from .sources.era5.constants import rg
import numpy as np


def _compute_ocean_only_mask(ds, z_max=5.0, lsm_max=0.2):
    """Computes a mask for cells that are considered to be over the
    ocean. The height of the surface is below z_max (in m, converted
    to geopotential using gravitational constant), and the lsm fraction
    is below lsm_max"""
    da_mask = (ds.z < z_max * rg) * (ds.lsm < lsm_max)
    da_mask.attrs["long_name"] = (
        "Ocean only (topography altitude < {z_max} [{ds.z.units}]"
        " and land-sea mask < {lsm_max} [{ds.lsm.units}])"
    )
    return da_mask


def _compute_all_mask(ds):
    da_mask = ds.z < ds.z + 1
    da_mask.attrs["long_name"] = "All data mask"
    return da_mask


def calc_mask(ds, mask_type):
    if mask_type == "ocean_only":
        return _compute_ocean_only_mask(ds=ds)
    elif mask_type == "all":
        return _compute_all_mask(ds=ds)
    else:
        raise NotImplementedError(mask_type)
