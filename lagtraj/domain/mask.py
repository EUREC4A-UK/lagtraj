from .sources.era5.constants import rg


def _compute_ocean_only_mask(ds, z_max=5.0, lsm_max=0.2):
    """Computes a mask for cells that are considered to be over the
    ocean. The height of the surface is below z_max (in m, converted to
    geopotential using gravitational constant), and the lsm fraction is
    below lsm_max. Here, lsm is the land-sea mask as a fraction (0 is
    all water, 1 is all land)."""
    da_mask = (ds.z < z_max * rg) * (ds.lsm < lsm_max)
    da_mask.attrs["long_name"] = (
        f"Ocean only (topography altitude < {z_max} [{ds.z.units}]"
        f" and land-sea mask < {lsm_max} [{ds.lsm.units}])"
    )
    return da_mask


def _compute_all_mask(ds):
    """Computes a mask that is true everywhere"""
    da_mask = ds.z < ds.z + 1
    da_mask.attrs["long_name"] = "All data mask"
    return da_mask


def _compute_land_only_mask(ds, lsm_min=0.8):
    """Computes a mask for cells that are considered to be over the
    land. The lsm fraction is above lsm_min. Here, lsm is the land-sea
    mask as a fraction (0 is all water, 1 is all land)."""
    da_mask = ds.lsm > lsm_min
    da_mask.attrs[
        "long_name"
    ] = f"Land only (land-sea mask > {lsm_min} [{ds.lsm.units}])"
    return da_mask


def calc_mask(ds, mask_type):
    if mask_type == "ocean_only":
        return _compute_ocean_only_mask(ds=ds)
    elif mask_type == "land_only":
        return _compute_land_only_mask(ds=ds)
    elif mask_type == "all":
        return _compute_all_mask(ds=ds)
    else:
        raise NotImplementedError(mask_type)
