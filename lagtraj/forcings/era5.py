from collections import namedtuple


ForcingSamplingDefinition = namedtuple(
    "ForcingSamplingDefinition",
    ["profile_method", "gradient_method", "averaging_width", "time_sampling_method",],
)


class InvalidLevelsDefinition(Exception):
    pass


def calculate_timestep(da_pt, ds_domain, sampling_method):
    """
    given the domain data in `ds_domain` calculate the forcing profile at a
    given point `da_pt` (da_pt is expected to have time, lat and lon positions)
    """
    ds_domain_timestep = ds_domain.sel(time=da_pt.time)
    if ds_domain_timestep.time.count() != 1:
        raise NotImplementedError(
            "Forcings based on era5 data cannot be"
            "interpolated between model data timesteps, and request time is"
            "outside"
        )

    if not da_pt.levels.units == "m":
        raise InvalidLevelsDefinition(
            "ERA5 data can currently only be generated on height levels"
        )

    if sampling_method.profile_method == "nearest_column":
        dvdt = ds_domain.differentiate(coord='time').sel(time=da_pt.time)
        ds_timestep_forcing = dvdt.sel(lat=da_pt.lat, lon=da_pt.lon, method="nearest")
    else:
        raise NotImplementedError(sampling_method)

    return ds_timestep_forcing
