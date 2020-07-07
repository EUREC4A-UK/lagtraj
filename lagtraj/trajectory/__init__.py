import isodate
from collections import namedtuple

from .. import build_data_path as build_data_path_global


TrajectoryOrigin = namedtuple("TrajectoryOrigin", ["lat", "lon", "datetime"])

TrajectoryDuration = namedtuple("TrajectoryDuration", ["forward", "backward"])

TrajectoryDefinition = namedtuple(
    "TrajectoryDefinition", ["domain", "duration", "origin", "name", "type", "timestep"]
)


INPUT_REQUIRED_FIELDS = {
    "trajectory_type": ["linear", "stationary"],
    "domain": str,
    "lat_origin": float,
    "lon_origin": float,
    "datetime_origin": isodate.parse_datetime,
    "forward_duration|backward_duration": isodate.parse_duration,
    "timestep": ("domain_data", isodate.parse_duration),
    "u_vel": [None, float],  # TODO: remove when velocity defn is improved
    "v_vel": [None, float],  # TODO: remove when velocity defn is imporved
}


def build_data_path(root_data_path, trajectory_name):
    data_path = build_data_path_global(
        root_data_path=root_data_path, data_type="trajectory"
    )
    return data_path / "{}.nc".format(trajectory_name)
