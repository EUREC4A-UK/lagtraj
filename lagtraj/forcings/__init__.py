from collections import namedtuple
from ..utils.interpolation.levels import ForcingLevelsDefinition  # noqa
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
    "ForcingDefinition", ["trajectory", "domain", "sampling", "levels", "name"],
)


def build_forcing_data_path(root_data_path, forcing_name, target_name=None):
    if target_name is None:
        filename = f"{forcing_name}.nc"
    else:
        filename = f"{forcing_name}.{target_name}.nc"
    return build_data_path(root_data_path=root_data_path, data_type="forcing") / (
        filename
    )
