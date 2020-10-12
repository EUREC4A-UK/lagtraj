from ..input_definitions import load
from . import (
    INPUT_REQUIRED_FIELDS,
    ConversionDefinition,
    ConversionParametersDefinition,
    ConversionMetadataDefinition,
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
        nudging_time_scalar_traj=conversion_params.get(
            "nudging_time_scalar_traj", None
        ),
        nudging_method_scalar_traj=conversion_params.get(
            "nudging_method_scalar_traj", None
        ),
        nudging_parameters_scalar_traj=conversion_params.get(
            "nudging_parameters_scalar_traj", None
        ),
        nudging_time_momentum_traj=conversion_params.get(
            "nudging_time_momentum_traj", None
        ),
        nudging_method_momentum_traj=conversion_params.get(
            "nudging_method_momentum_traj", None
        ),
        nudging_parameters_momentum_traj=conversion_params.get(
            "nudging_parameters_momentum_traj", None
        ),
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
