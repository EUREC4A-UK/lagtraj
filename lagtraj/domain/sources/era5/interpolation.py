import numbers

import xarray as xr
import numpy as np
from ....utils.interpolation import steffen_3d
from .utils import calculate_heights_and_pressures


def interpolate_to_height_levels(ds_model_levels, height, mask_method="sea"):
    """Converts ERA5 model level data to data on height levels
    using Steffen interpolation"""
    height = np.atleast_1d(height)

    if "height_f" not in ds_model_levels:
        ds_extra = calculate_heights_and_pressures(ds_model_levels)
        ds_model_levels = ds_model_levels.merge(ds_extra)

    if isinstance(height[0], numbers.Integral):
        raise Exception("Heights need to be floating numbers, rather than integers")

    # `steffen_3d` expects the data to have the shape
    # (level, lat, lon) so we explicitly expand the dataset here and squeeze
    # later
    ds_ = ds_model_levels
    required_dims = ("time", "level", "lat", "lon")
    missing_dims = list(filter(lambda d: d not in ds_.dims, required_dims))
    ds_ = ds_.expand_dims(missing_dims).transpose(*required_dims)

    datasets_on_height_levels = []
    for t in ds_.time.values:
        ds_timestep_model_levels = ds_.sel(time=t)

        # construct coords for output dataset (on height levels)
        coords = dict(ds_timestep_model_levels.coords)
        del coords["level"]
        coords["height"] = xr.DataArray(
            height, attrs={"long_name": "altitude", "units": "metres"}
        )
        ds_timestep_height_levels = xr.Dataset(coords=coords)

        # Points over sea that have height above zero but below 5m
        # Here, profiles are extended to 0m
        z_min_2d = ds_timestep_model_levels.height_h.min(dim="level")
        lsm_2d = ds_timestep_model_levels.lsm
        sea_mask = (z_min_2d < 5.0) * (z_min_2d > 1.0e-6) * (lsm_2d < 0.2)
        z_min_surface = (
            ds_timestep_model_levels.height_h.isel(level=-1)
            .where(~sea_mask, other=-1.0e-6)
            .values
        )
        z_max_surface = ds_timestep_model_levels.height_f.isel(level=0).values

        # pressure levels are monitonically in the opposite order, to ensure
        # the heights end up in increasing order we flip the input and height
        # level arrays
        def get_height_reversed_values(da):
            assert da.dims == ("level", "lat", "lon")
            return da[::-1, :, :].values

        for v in ds_timestep_model_levels.data_vars:
            da_v = ds_timestep_model_levels[v]
            if "level" not in da_v.dims:
                # doesn't need interpolating so we just copy it across
                da_v = ds_timestep_model_levels[v]
            else:
                field_p_levels = ds_timestep_model_levels[v]
                interp_kwargs = dict(
                    v_in=get_height_reversed_values(field_p_levels),
                    z_out=height,
                    z_min_surface=z_min_surface,
                    z_max_surface=z_max_surface,
                )

                if v in ["height_h", "p_h"]:
                    interp_kwargs["z_in"] = get_height_reversed_values(
                        ds_timestep_model_levels.height_h
                    )
                else:
                    interp_kwargs["z_in"] = get_height_reversed_values(
                        ds_timestep_model_levels.height_f
                    )

                if v in ["p_h", "p_f", "height_h", "height_f"]:
                    # height and pressure values should be extrapolated using
                    # the gradient at the lower and upper boundary
                    interp_kwargs["lower_extrapolation_with_gradient"] = True
                    interp_kwargs["upper_extrapolation_with_gradient"] = True

                da_v_height_levels = xr.DataArray(
                    steffen_3d(**interp_kwargs), dims=("height", "lat", "lon")
                )
                ds_timestep_height_levels[v] = da_v_height_levels

        datasets_on_height_levels.append(ds_timestep_height_levels)

    return xr.concat(datasets_on_height_levels, dim="time")
