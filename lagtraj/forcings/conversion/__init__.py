import numpy as np
from collections import namedtuple


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
    nudging_u=(0, np.nan, float),
    nudging_v=(0, np.nan, float),
    nudging_temp=(0, np.nan, float),
    nudging_theta=(0, np.nan, float),
    nudging_thetal=(0, np.nan, float),
    nudging_qv=(0, np.nan, float),
    nudging_qt=(0, np.nan, float),
    nudging_rv=(0, np.nan, float),
    nudging_rt=(0, np.nan, float),
    surfaceType=["ocean", "land", "mixed"],
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
    inversion_nudging=(None, int),
    inversion_nudging_height_above=(None, float),
    inversion_nudging_transition=(None,float),
    inversion_nudging_time=(None,float)
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
        "inversion_nudging",
        "inversion_nudging_height_above",
        "inversion_nudging_transition",
        "inversion_nudging_time",
        "surfaceType",
        "surfaceForcing",
        "surfaceForcingWind",
        "nudging_parameters_scalar_traj",
        "nudging_parameters_momentum_traj",
    ],
)

ConversionNudgingDefinition = namedtuple(
    "ConversionNudgingDefinition", ["method", "time", "height", "transition",],
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

from . import targets  # noqa
from ...utils.interpolation.levels import (
    LevelsDefinition as ConversionLevelsDefinition,
)  # noqa
from .process import export_for_target
