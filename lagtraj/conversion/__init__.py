from collections import namedtuple
from .utils.levels import ConversionLevelsDefinition
from .. import build_data_path


INPUT_REQUIRED_FIELDS = dict(
    export_format=str,
    levels_method=(None, str),
    levels_number=(None, int),
    levels_dzmin=(None, float),
    levels_ztop=(None, float),
    nudging_time=float,
)

ConversionDefinition = namedtuple(
    "ConversionDefinition", ["export_format", "levels", "name", "nudging_time"],
)


def build_conversion_data_path(root_data_path, forcing_name, conversion_name):
    return build_data_path(root_data_path=root_data_path, data_type="conversion") / (
        forcing_name + "_" + conversion_name + ".nc"
    )
