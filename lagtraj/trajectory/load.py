import datetime

import xarray as xr

from ..input_definitions import load
from . import INPUT_REQUIRED_FIELDS
from . import (
    TrajectoryDuration,
    TrajectoryOrigin,
    build_data_path,
    TrajectoryDefinition,
)


def load_definition(root_data_path, name):
    params = load.load_definition(
        root_data_path=root_data_path,
        input_name=name,
        input_type="trajectory",
        required_fields=INPUT_REQUIRED_FIELDS,
    )

    origin = TrajectoryOrigin(
        lat=params["lat_origin"],
        lon=params["lon_origin"],
        datetime=params["datetime_origin"],
    )

    duration = TrajectoryDuration(
        forward=params.get("forward_duration", datetime.timedelta()),
        backward=params.get("backward_duration", datetime.timedelta()),
    )

    return TrajectoryDefinition(
        domain=params["domain"],
        origin=origin,
        duration=duration,
        type=params["trajectory_type"],
        name=params["name"],
        timestep=params["timestep"],
    )


def load_data(root_data_path, name):
    trajectory_data_path = build_data_path(
        root_data_path=root_data_path, trajectory_name=name
    )
    ds = xr.open_dataset(trajectory_data_path)
    return ds
