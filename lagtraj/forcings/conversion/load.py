import shutil


from ...input_definitions import (
    load, build_input_definition_path, examples as input_examples
)
from . import (
    INPUT_REQUIRED_FIELDS,
    ConversionDefinition,
    ConversionParametersDefinition,
    ConversionMetadataDefinition,
    ConversionLevelsDefinition,
    ConversionNudgingDefinition,
)


def load_definition(root_data_path, forcing_name, conversion_name, target_name):
    conversion_defn_path = build_input_definition_path(
        root_data_path=root_data_path,
        input_name=forcing_name,
        input_type="forcing",
        input_subtype=target_name,
    )

    # first we look if there's a conversion defined specifically for
    # targeting a forcing with name `forcing_name` to the model with name
    # `target_name` (e.g. kpt or dephy).  If not we try to see if there's a
    # default conversion-file we can make a copy of and use
    if not conversion_defn_path.exists():
        default_conversion_defn_path = build_input_definition_path(
            root_data_path=root_data_path,
            input_name=target_name,
            input_type="forcing_conversions",
        )
        if not default_conversion_defn_path.exists():
            available_default_conversions = input_examples.get_available(
                input_types="forcing_conversions"
            )
            s_avail = ", ".join(available_default_conversions.keys())
            raise Exception(
                f"Couldn't find a forcing conversion file in `{conversion_defn_path}`"
                f" for converting the `{forcing_name}` forcing to target the"
                f" `{target_name}` model, and a default input-definition for"
                f" how to target the `{target_name}` isn't included currently"
                f" with lagtraj. Please use one of the currently bundled"
                f" conversions included with lagtraj ({s_avail})"
                f" by targeting one of those models."
            )
        else:
            shutil.copy(default_conversion_defn_path, conversion_defn_path)
            print(
                f"Parameters for how to target the `{target_name}` model for the "
                f"`{forcing_name}` forcing weren't found in `{conversion_defn_path}` "
                f"and so default parameters for `{target_name}` were copied to "
                f"{conversion_defn_path}. Please change these parameter as needed "
                "and rerun the conversion."
            )
    conversion_params = load.load_definition(
        root_data_path=root_data_path,
        input_name=forcing_name,
        input_type="forcing",
        input_subtype=target_name,
        required_fields=INPUT_REQUIRED_FIELDS,
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
