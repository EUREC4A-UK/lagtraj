from pathlib import Path
import os


__version__ = "0.1.0"


# Optional numba dependency
try:
    from numba import njit

    print("Running with numba")
except ImportError:

    def njit(numba_function):
        """Dummy numba function"""
        return numba_function

    print("Running without numba")


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
