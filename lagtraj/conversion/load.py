from ..input_definitions import load
from . import (
    INPUT_REQUIRED_FIELDS,
    ConversionDefinition,
    ConversionLevelsDefinition,
)


def load_definition(root_data_path, conversion_name):
    conversion_params = load.load_definition(
        root_data_path=root_data_path,
        input_name=conversion_name,
        input_type="conversion",
        required_fields=INPUT_REQUIRED_FIELDS,
    )

    conversion_levels_definition = ConversionLevelsDefinition(
        n_levels=conversion_params.get("levels_number", None),
        z_top=conversion_params.get("levels_ztop", None),
        dz_min=conversion_params.get("levels_dzmin", None),
        method=conversion_params.get("levels_method", None),
    )

    conversion_definition = ConversionDefinition(
        export_format=conversion_params["export_format"],
        nudging_time=conversion_params["nudging_time"],
        levels=conversion_levels_definition,
        name=conversion_params["name"],
    )

    return conversion_definition
