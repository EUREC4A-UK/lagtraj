from .constants import rg


def era5_mask(ds_to_mask, dictionary):
    """Returns a lat-lon mask"""
    # Only use ocean points, ensure it can be used after before or after array extensions
    if dictionary["mask"] == "ocean":
        mask = (ds_to_mask["z"] < 5.0 * rg) * (ds_to_mask["lsm"].values < 0.2)
    return mask
