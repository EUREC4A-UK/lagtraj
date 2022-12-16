from collections import namedtuple

from ...utils.interpolation.levels import (  # noqa
    LevelsDefinition as ConversionLevelsDefinition,
)

INPUT_REQUIRED_FIELDS = dict(
    export_format=str,
    levels_method=(None, str),
    levels_number=(None, int),
    levels_dzmin=(None, float),
    levels_ztop=(None, float),
    comment=(None, str),
    campaign=(None, str),
    source_domain=(None, str),
    reference=(None, str),
    # AUTHOR=CREATOR
    author=(None, str),
    modifications=(None, str),
    # CASE=FLIGHT
    case=(None, str),
    adv_temp=[0, 1],
    adv_theta=[0, 1],
    adv_thetal=[0, 1],
    adv_qv=[0, 1],
    adv_qt=[0, 1],
    adv_rv=[0, 1],
    adv_rt=[0, 1],
    rad_temp=[0, 1, "adv"],
    rad_theta=[0, 1, "adv"],
    rad_thetal=[0, 1, "adv"],
    forc_omega=[0, 1],
    forc_w=[0, 1],
    forc_geo=[0, 1],
    # Allow integers only to conform to DEPHY standard, for more complex
    # options, use the nudging "methods" below
    nudging_u=(int),
    nudging_v=(int),
    nudging_temp=(int),
    nudging_theta=(int),
    nudging_thetal=(int),
    nudging_qv=(int),
    nudging_qt=(int),
    nudging_rv=(int),
    nudging_rt=(int),
    surfaceType=["ocean", "land", "mixed", "landice"],
    surfaceForcing=["ts", "Flux", "surfaceFlux"],
    surfaceForcingWind=["z0", "ustar", "z0_traj"],
    nudging_method_scalar_traj=(None, str),
    nudging_time_scalar_traj=(None, float),
    nudging_height_scalar_traj=(None, float),
    nudging_transition_scalar_traj=(None, float),
    nudging_method_momentum_traj=(None, str),
    nudging_time_momentum_traj=(None, float),
    nudging_height_momentum_traj=(None, float),
    nudging_transition_momentum_traj=(None, float),
    # Nudging above inversion, for EUREC4A cases
    inversion_nudging=(None, int),
    inversion_nudging_height_above=(None, float),
    inversion_nudging_transition=(None, float),
    inversion_nudging_time=(None, float),
    # Remove tendencies and set geowind to wind at high levels
    wind_at_high_levels_correction=(None, int),
    wind_at_high_levels_correction_pressure_above=(None, float),
    wind_at_high_levels_correction_transition=(None, float),
)

ConversionDefinition = namedtuple(
    "ConversionDefinition",
    ["export_format", "levels", "name", "metadata", "parameters"],
)

ConversionParametersDefinition = namedtuple(
    "ConversionParametersDefinition",
    [
        "adv_temp",
        "adv_theta",
        "adv_thetal",
        "adv_qv",
        "adv_qt",
        "adv_rv",
        "adv_rt",
        "rad_temp",
        "rad_theta",
        "rad_thetal",
        "forc_omega",
        "forc_w",
        "forc_geo",
        "nudging_u",
        "nudging_v",
        "nudging_temp",
        "nudging_theta",
        "nudging_thetal",
        "nudging_qv",
        "nudging_qt",
        "nudging_rv",
        "nudging_rt",
        "surfaceType",
        "surfaceForcing",
        "surfaceForcingWind",
        "nudging_parameters_scalar_traj",
        "nudging_parameters_momentum_traj",
        "inversion_nudging",
        "inversion_nudging_height_above",
        "inversion_nudging_transition",
        "inversion_nudging_time",
        "wind_at_high_levels_correction",
        "wind_at_high_levels_correction_pressure_above",
        "wind_at_high_levels_correction_transition",
    ],
)

ConversionNudgingDefinition = namedtuple(
    "ConversionNudgingDefinition",
    [
        "method",
        "time",
        "height",
        "transition",
    ],
)

ConversionMetadataDefinition = namedtuple(
    "ConversionMetadataDefinition",
    [
        "comment",
        "campaign",
        "source_domain",
        "reference",
        "author",
        "modifications",
        "case",
    ],
)
