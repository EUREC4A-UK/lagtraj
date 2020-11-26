from collections import namedtuple
from .utils.levels import ForcingLevelsDefinition  # noqa
from .profile_calculation import ForcingSamplingDefinition  # noqa
from .. import build_data_path


INPUT_REQUIRED_FIELDS = dict(
    trajectory=str,
    domain=str,
    gradient_method=str,
    averaging_width=float,
    levels_method=str,
    levels_number=int,
    levels_dzmin=(None, float),
    levels_ztop=float,
    time_sampling_method=(None, "str"),
    sampling_mask=[None, "ocean_only"],
)


ForcingDefinition = namedtuple(
    "ForcingDefinition",
    ["trajectory", "domain", "sampling", "levels", "name", "version"],
)


def build_forcing_data_path(root_data_path, forcing_name):
    return build_data_path(root_data_path=root_data_path, data_type="forcing") / (
        forcing_name + ".nc"
    )
