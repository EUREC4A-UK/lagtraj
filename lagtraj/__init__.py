from pathlib import Path
import os

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


DATA_TYPE_PLURAL = dict(trajectory="trajectories", domain="domains", forcing="forcings")


def build_data_path(root_data_path, data_type):
    data_type_plural = DATA_TYPE_PLURAL[data_type]
    return Path(root_data_path) / data_type_plural

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
