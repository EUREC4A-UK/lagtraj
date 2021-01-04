import isodate
from collections import namedtuple

from .. import build_data_path as build_data_path_global


TrajectoryOrigin = namedtuple("TrajectoryOrigin", ["lat", "lon", "datetime"])

TrajectoryDuration = namedtuple("TrajectoryDuration", ["forward", "backward"])

TrajectoryDefinition = namedtuple(
    "TrajectoryDefinition",
    ["domain", "duration", "origin", "name", "type", "timestep", "extra_kwargs"],
)


INPUT_REQUIRED_FIELDS = {
    "trajectory_type": ["linear", "eulerian", "lagrangian"],
    # domain should only be given when creating a lagrangian trajectory
    "domain": dict(requires=dict(trajectory_type="lagrangian"), choices=str),
    "lat_origin": float,
    "lon_origin": float,
    "datetime_origin": isodate.parse_datetime,
    "forward_duration|backward_duration": isodate.parse_duration,
    "timestep": ("domain_data", isodate.parse_duration),
    # only linear trajectories need to have their velocity prescribed
    "u_vel": dict(requires=dict(trajectory_type="linear"), choices=float),
    "v_vel": dict(requires=dict(trajectory_type="linear"), choices=float),
    # velocity method is only relevant when making a lagrangian trajectories
    "velocity_method": dict(
        requires=dict(trajectory_type="lagrangian"),
        choices=[
            "single_height_level",
            "single_pressure_level",
            "lower_troposphere_humidity_weighted"
        ],
    ),
    "velocity_method_height": dict(
        requires=dict(velocity_method="single_height_level"),
        choices=float,
    ),
    "velocity_method_pressure": dict(
        requires=dict(velocity_method="single_pressure_level"),
        choices=float,
    ),
}


def build_data_path(root_data_path, trajectory_name):
    data_path = build_data_path_global(
        root_data_path=root_data_path, data_type="trajectory"
    )
    return data_path / "{}.nc".format(trajectory_name)
