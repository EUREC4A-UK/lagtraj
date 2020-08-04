from ..input_definitions import load
from . import (
    INPUT_REQUIRED_FIELDS,
    ForcingDefinition,
    ForcingLevelsDefinition,
    ForcingSamplingDefinition,
)


def load_definition(root_data_path, forcing_name):
    forcing_params = load.load_definition(
        root_data_path=root_data_path,
        input_name=forcing_name,
        input_type="forcing",
        required_fields=INPUT_REQUIRED_FIELDS,
    )

    forcing_levels_definition = ForcingLevelsDefinition(
        n_levels=forcing_params["levels_number"],
        z_top=forcing_params["levels_ztop"],
        dz_min=forcing_params.get("levels_dzmin", None),
        method=forcing_params["levels_method"],
    )

    sampling_definition = ForcingSamplingDefinition(
        gradient_method=forcing_params["gradient_method"],
        averaging_width=forcing_params["averaging_width"],
        time_sampling_method=forcing_params.get("time_sampling_method", "domain_data"),
    )

    forcing_definition = ForcingDefinition(
        trajectory=forcing_params["trajectory"],
        domain=forcing_params["domain"],
        sampling=sampling_definition,
        levels=forcing_levels_definition,
        name=forcing_params["name"],
    )

    return forcing_definition
