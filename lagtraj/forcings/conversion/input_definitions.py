from collections import namedtuple

from ...utils.interpolation.levels import (  # noqa
    LevelsDefinition as ConversionLevelsDefinition,
)
from .input_definitions_nudging import NUDGING_REQUIRED_FIELDS

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
    surfaceType=["ocean", "land", "mixed"],
    surfaceForcing=["ts", "Flux", "surfaceFlux"],
    surfaceForcingWind=["z0", "ustar", "z0_traj"],
)

# add parameters for nudging of forcing profiles for scalars and momentum
INPUT_REQUIRED_FIELDS.update(NUDGING_REQUIRED_FIELDS)

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
        "surfaceType",
        "surfaceForcing",
        "surfaceForcingWind",
    ]
    + list(NUDGING_REQUIRED_FIELDS.keys()),
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
