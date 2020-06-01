import numpy as np
from scipy.optimize import ridder

import xarray as xr

# Utilities related to level spacing and interpolation
# Except steffen interpolation, which is in separate pyx file


def make_levels(forcings_dict):
    top = forcings_dict["levels_top"]
    n_levels = forcings_dict["levels_number"]
    if forcings_dict["levels_strategy"] == "linear":
        dz_min = forcings_dict["levels_dzmin"]
        levels = make_exponential_levels(dz_min, top, n_levels)
    elif forcings_dict["levels_strategy"] == "exponential":
        levels = make_linear_levels(top, n_levels)
    else:
        raise Exception("Unknown strategy for making levels")
    return xr.DataArray(
        levels, dims=("level",), attrs=dict(long_name="height", units="m")
    )


def exponential_levels_error(x, dz_min, top, n_levels):
    return top - dz_min * ((1.0 - x ** (n_levels - 1)) / (1.0 - x))


# Makes exponentially distributed levels
def make_exponential_levels(dz_min, top, n_levels):
    level_rate = ridder(
        exponential_levels_error, 1.000000000000001, 2.0, args=(dz_min, top, n_levels)
    )
    level_list = [0.0]
    for level in range(1, n_levels):
        level_list = level_list + [
            dz_min * ((1.0 - level_rate ** level) / (1.0 - level_rate))
        ]
    level_list[n_levels - 1] = top
    return np.array(level_list)


def make_linear_levels(top, n_levels):
    level_array = np.linspace(0.0, top, n_levels)
    return level_array
