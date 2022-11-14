import os
import warnings
from pathlib import Path

__version__ = "0.1.1"


try:
    from numba import njit
except ImportError:

    def njit(fn):
        """Dummy numba function"""
        return fn

    warnings.warn(
        "Running without numba. Computation will be severely slowed"
        "down. Please install `numba` to use it with `lagtraj`"
    )


# by default we store data relative to where lagtraj is invoked from
DEFAULT_ROOT_DATA_PATH = Path(os.getcwd()) / "data"


DATA_TYPE_PLURAL = dict(
    trajectory="trajectories",
    domain="domains",
    forcing="forcings",
    forcing_conversion="forcings_conversion",
)


def build_data_path(root_data_path, data_type):
    data_type_plural = DATA_TYPE_PLURAL[data_type]
    return Path(root_data_path) / data_type_plural
