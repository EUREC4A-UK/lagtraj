from collections import namedtuple

from .. import build_data_path
from .sources import (
    interpolate_to_height_levels,
    interpolate_to_pressure_levels,
)  # noqa
from .sources import calc_auxiliary_variable  # noqa

LatLonBoundingBox = namedtuple(
    "LatLonBoundingBox", ["lat_min", "lat_max", "lon_min", "lon_max"]
)

LatLonSamplingResolution = namedtuple("LatLonSamplingResolution", ["lat", "lon"])


INPUT_REQUIRED_FIELDS = dict(
    source=str,
    lat_min=float,
    lat_max=float,
    lon_min=float,
    lon_max=float,
    lat_samp=float,
    lon_samp=float,
)


def build_domain_data_path(root_data_path, domain_name):
    return build_data_path(root_data_path=root_data_path, data_type="domain") / (
        f"{domain_name}_data"
    )
