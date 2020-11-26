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

    extra_kwargs = {}
    if "u_vel" in params or "v_vel" in params:
        if not ("u_vel" in params and "v_vel" in params):
            raise Exception(
                "Both `u_vel` and `v_vel` should be defined when"
                " creating a linear trajectory"
            )
        extra_kwargs["U"] = (params["u_vel"], params["v_vel"])

    if "velocity_method" in params:
        extra_kwargs["velocity_method"] = params["velocity_method"]

    if "velocity_method_height" in params:
        extra_kwargs["velocity_method_kwargs"] = dict(
            height=params["velocity_method_height"]
        )

    if "velocity_method_pressure" in params:
        extra_kwargs["velocity_method_kwargs"] = dict(
            pressure=params["velocity_method_pressure"]
        )

    return TrajectoryDefinition(
        domain=params["domain"],
        origin=origin,
        duration=duration,
        type=params["trajectory_type"],
        name=params["name"],
        timestep=params["timestep"],
        extra_kwargs=extra_kwargs,
        version=params.get("version", "unversioned"),
    )


def load_data(root_data_path, name):
    trajectory_definition = load_definition(root_data_path=root_data_path, name=name)

    trajectory_data_path = build_data_path(
        root_data_path=root_data_path, trajectory_name=name
    )
    ds = xr.open_dataset(trajectory_data_path)

    if not "version" in ds.attrs:
        raise Exception(
            f"The trajectory stored in `{trajectory_data_path}` "
            "doesn't `version` attribute set. Please delete and "
            "recreate the trajectory or set the correct version. "
            "Was expecting to find trajectory with version "
            f"`{trajectory_definition['version']}`"
        )

    if ds.version != trajectory_definition["version"]:
        raise Exception(
            f"The version of the trajectory stored in `{trajectory_data_path}` "
            "doesn't match the version in the input definition yaml file "
            f"for a trajectory named `{name}`. Please delete the trajectory "
            " and recreate it"
        )

    return ds
