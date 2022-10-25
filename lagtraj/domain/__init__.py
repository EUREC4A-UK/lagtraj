from collections import namedtuple

from .. import build_data_path
from ..input_definitions.examples import LAGTRAJ_EXAMPLES_PATH_PREFIX
from .sources import calc_auxiliary_variable  # noqa
from .sources import (  # noqa
    interpolate_to_height_levels,
    interpolate_to_pressure_levels,
)

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
    # we need to strip the `lagtraj://` prefix before we construct the path
    # since the data is stored locally
    if domain_name.startswith(LAGTRAJ_EXAMPLES_PATH_PREFIX):
        domain_name = domain_name[len(LAGTRAJ_EXAMPLES_PATH_PREFIX) :]

    return build_data_path(root_data_path=root_data_path, data_type="domain") / (
        f"{domain_name}_data"
    )
