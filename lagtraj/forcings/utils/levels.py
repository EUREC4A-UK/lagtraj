"""
Utilities related to level spacing and interpolation
Except steffen interpolation, which is in separate pyx file
"""
from collections import namedtuple

import numpy as np
from scipy.optimize import ridder

import xarray as xr


ForcingLevelsDefinition = namedtuple(
    "ForcingLevelsDefinition", ["method", "n_levels", "z_top", "dz_min"]
)


def make_levels(method, n_levels, z_top, dz_min=None):
    if method == "linear":
        levels = _make_exponential_levels(dz_min, z_top, n_levels)
    elif method == "exponential":
        levels = _make_linear_levels(z_top, n_levels)
    else:
        raise Exception("Unknown strategy for making levels")
    return xr.DataArray(
        levels,
        dims=("level",),
        attrs=dict(long_name="height", units="m", method="linear"),
    )


def exponential_levels_error(x, dz_min, z_top, n_levels):
    return z_top - dz_min * ((1.0 - x ** (n_levels - 1)) / (1.0 - x))


def _make_exponential_levels(dz_min, z_top, n_levels):
    """
    Makes exponentially distributed levels
    """
    if dz_min is None:
        raise Exception("For exponential levels `dz_min` must be a float")
    level_rate = ridder(
        exponential_levels_error, 1.000000000000001, 2.0, args=(dz_min, z_top, n_levels)
    )
    level_list = [0.0]
    for level in range(1, n_levels):
        level_list = level_list + [
            dz_min * ((1.0 - level_rate ** level) / (1.0 - level_rate))
        ]
    level_list[n_levels - 1] = z_top
    return np.array(level_list)


def _make_linear_levels(top, n_levels):
    level_array = np.linspace(0.0, top, n_levels)
    return level_array
