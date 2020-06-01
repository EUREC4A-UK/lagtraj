from collections import namedtuple

from .. import build_data_path

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
        domain_name + "_data"
    )
