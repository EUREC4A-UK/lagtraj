from ...input_definitions import (
    load,
    build_input_definition_path,
)
from ...input_definitions.examples import LAGTRAJ_EXAMPLES_PATH_PREFIX
from . import (
    INPUT_REQUIRED_FIELDS,
    ConversionDefinition,
    ConversionParametersDefinition,
    ConversionMetadataDefinition,
    ConversionLevelsDefinition,
    ConversionNudgingDefinition,
)


def _get_definition_parameters(root_data_path, forcing_name, conversion_name):
    # if the user has requested a conversion target with the `lagtraj://`
    # prefix then they are asking for the default parameters for targeting
    # the `conversion_name` model. We will try to load that and then if
    # successful copy this input definition (for the conversion
    # specifically) to their local path (in the folder with the forcing
    # input definition)
    if conversion_name.startswith(LAGTRAJ_EXAMPLES_PATH_PREFIX):
        params_defn_expected_local_path = build_input_definition_path(
            root_data_path=root_data_path,
            input_name=forcing_name,
            input_type="forcing",
            input_subtype=conversion_name.replace(LAGTRAJ_EXAMPLES_PATH_PREFIX, ""),
        )
        conversion_defn = load.load_definition(
            root_data_path=root_data_path,
            input_name=conversion_name,
            input_type="forcing_conversion",
            required_fields=INPUT_REQUIRED_FIELDS,
            expected_local_path=params_defn_expected_local_path,
        )
        return conversion_defn
    else:
        return load.load_definition(
            root_data_path=root_data_path,
            input_name=forcing_name,
            input_type="forcing",
            input_subtype=conversion_name,
            required_fields=INPUT_REQUIRED_FIELDS,
        )


def load_definition(root_data_path, forcing_name, conversion_name):
    conversion_params = _get_definition_parameters(
        root_data_path=root_data_path,
        forcing_name=forcing_name,
        conversion_name=conversion_name,
    )

    conversion_levels_definition = ConversionLevelsDefinition(
        n_levels=conversion_params.get("levels_number", None),
        z_top=conversion_params.get("levels_ztop", None),
        dz_min=conversion_params.get("levels_dzmin", None),
        method=conversion_params.get("levels_method", None),
    )

    conversion_nudging_momentum_definition = ConversionNudgingDefinition(
        method=conversion_params.get("nudging_method_momentum_traj", None),
        time=conversion_params.get("nudging_time_momentum_traj", None),
        height=conversion_params.get("nudging_height_momentum_traj", None),
        transition=conversion_params.get("nudging_transition_momentum_traj", None),
    )

    conversion_nudging_scalar_definition = ConversionNudgingDefinition(
        method=conversion_params.get("nudging_method_scalar_traj", None),
        time=conversion_params.get("nudging_time_scalar_traj", None),
        height=conversion_params.get("nudging_height_scalar_traj", None),
        transition=conversion_params.get("nudging_transition_scalar_traj", None),
    )

    conversion_parameters_definition = ConversionParametersDefinition(
        adv_temp=conversion_params["adv_temp"],
        adv_theta=conversion_params["adv_theta"],
        adv_thetal=conversion_params["adv_thetal"],
        adv_qv=conversion_params["adv_qv"],
        adv_qt=conversion_params["adv_qt"],
        adv_rv=conversion_params["adv_rv"],
        adv_rt=conversion_params["adv_rt"],
        rad_temp=conversion_params["rad_temp"],
        rad_theta=conversion_params["rad_theta"],
        rad_thetal=conversion_params["rad_thetal"],
        forc_omega=conversion_params["forc_omega"],
        forc_w=conversion_params["forc_w"],
        forc_geo=conversion_params["forc_geo"],
        nudging_u=conversion_params["nudging_u"],
        nudging_v=conversion_params["nudging_v"],
        nudging_temp=conversion_params["nudging_temp"],
        nudging_theta=conversion_params["nudging_theta"],
        nudging_thetal=conversion_params["nudging_thetal"],
        nudging_qv=conversion_params["nudging_qv"],
        nudging_qt=conversion_params["nudging_qt"],
        nudging_rv=conversion_params["nudging_rv"],
        nudging_rt=conversion_params["nudging_rt"],
        surfaceType=conversion_params["surfaceType"],
        surfaceForcing=conversion_params["surfaceForcing"],
        surfaceForcingWind=conversion_params["surfaceForcingWind"],
        nudging_parameters_momentum_traj=conversion_nudging_momentum_definition,
        nudging_parameters_scalar_traj=conversion_nudging_scalar_definition,
        inversion_nudging=conversion_params.get("inversion_nudging", None),
        inversion_nudging_height_above=conversion_params.get("inversion_nudging_height_above", None),
        inversion_nudging_transition=conversion_params.get("inversion_nudging_transition", None),
        inversion_nudging_time=conversion_params.get("inversion_nudging_time", None),
    )

    conversion_metadata_definition = ConversionMetadataDefinition(
        comment=conversion_params.get("comment", None),
        campaign=conversion_params.get("campaign", None),
        source_domain=conversion_params.get("source_domain", None),
        reference=conversion_params.get("reference", None),
        author=conversion_params.get("author", None),
        modifications=conversion_params.get("modifications", None),
        case=conversion_params.get("case", None),
    )

    conversion_definition = ConversionDefinition(
        export_format=conversion_params["export_format"],
        levels=conversion_levels_definition,
        parameters=conversion_parameters_definition,
        metadata=conversion_metadata_definition,
        name=conversion_params["name"],
    )

    return conversion_definition
