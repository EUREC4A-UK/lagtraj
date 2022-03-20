from collections import namedtuple

from .. import build_data_path
from ..input_definitions.examples import LAGTRAJ_EXAMPLES_PATH_PREFIX
from ..utils.interpolation.levels import (  # noqa
    LevelsDefinition as ForcingLevelsDefinition,
)
from .profile_calculation import ForcingSamplingDefinition  # noqa

INPUT_REQUIRED_FIELDS = dict(
    trajectory=str,
    domain=str,
    gradient_method=str,
    advection_velocity_sampling_method=str,
    averaging_width=float,
    levels_method=str,
    levels_number=int,
    levels_dzmin=(None, float),
    levels_ztop=float,
    time_sampling_method=(None, "str"),
    sampling_mask=["all", "land_only", "ocean_only"],
)


ForcingDefinition = namedtuple(
    "ForcingDefinition",
    ["trajectory", "domain", "sampling", "levels", "name", "version"],
)


def build_forcing_data_path(root_data_path, forcing_name, conversion_name=None):
    """
    Create the output path for a forcing with name `forcing_name` in the path
    `root_data_path`, optionally targetting a model by name `conversion_name`
    """
    # we need to strip the `lagtraj://` prefix before we construct the path
    # since the data is stored locally
    if forcing_name.startswith(LAGTRAJ_EXAMPLES_PATH_PREFIX):
        forcing_name = forcing_name[len(LAGTRAJ_EXAMPLES_PATH_PREFIX) :]

    if conversion_name is not None and conversion_name.startswith(
        LAGTRAJ_EXAMPLES_PATH_PREFIX
    ):
        conversion_name = conversion_name[len(LAGTRAJ_EXAMPLES_PATH_PREFIX) :]

    if conversion_name is None:
        filename = f"{forcing_name}.nc"
    else:
        filename = f"{forcing_name}.{conversion_name}.nc"

    return build_data_path(root_data_path=root_data_path, data_type="forcing") / (
        filename
    )
