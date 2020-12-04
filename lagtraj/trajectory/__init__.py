import isodate
from collections import namedtuple

from .. import build_data_path as build_data_path_global
from ..input_definitions import InvalidInputDefinition


TrajectoryOrigin = namedtuple("TrajectoryOrigin", ["lat", "lon", "datetime"])

TrajectoryDuration = namedtuple("TrajectoryDuration", ["forward", "backward"])

TrajectoryDefinition = namedtuple(
    "TrajectoryDefinition",
    [
        "domain",
        "duration",
        "origin",
        "name",
        "type",
        "timestep",
        "extra_kwargs",
        "version",
    ],
)


def _validate_domain_param(param_name, input_params):
    """
    domain should only be given for `eulerian` and `linear` trajectories if the
    timesteps are to be extracted from it
    """
    if input_params["trajectory_type"] == "lagrangian":
        if param_name not in input_params:
            raise InvalidInputDefinition(
                "Lagrangian trajectories require a domain for the integration. "
                "Please set the `domain` variable to name the source data domain."
            )
        elif not isinstance(input_params.get(param_name, None), str):
            raise InvalidInputDefinition(
                "The parameter `domain` should be a string defining the named domain "
                "to use for the trajectory calculation."
            )
        return input_params[param_name]

    elif input_params["trajectory_type"] in ["eulerian", "linear"]:
        if (
            param_name in input_params
            and input_params.get("timestep", None) != "domain_data"
        ):
            raise InvalidInputDefinition(
                "The `domain` parameter should only be defined if the trajectory "
                "timesteps are to be taken from it (`timestep`==`domain_data`)"
            )
        return None

    else:
        raise NotImplementedError(input_params["trajectory_type"])


def _validate_velocity_param(param_name, input_params):
    if input_params["trajectory_type"] == "linear":
        if param_name not in input_params:
            raise InvalidInputDefinition(
                f"To create a linear trajectory you must define `{param_name}`"
            )
        else:
            param_value = input_params[param_name]
            if not (isinstance(param_value, float) or isinstance(param_value, int)):
                raise InvalidInputDefinition(
                    f"The trajectory velocity `{param_name}` should be a number "
                    f"(it is currently `{input_params[param_name]}`)"
                )
            else:
                return param_value
    elif input_params["trajectory_type"] in ["eulerian", "lagrangian"]:
        if param_name in input_params:
            raise InvalidInputDefinition(
                f"The parameter `{param_name}` shouldn't be defined for "
                f"`{input_params['trajectory_type']}` trajectories"
            )


INPUT_REQUIRED_FIELDS = {
    "trajectory_type": ["linear", "eulerian", "lagrangian"],
    "domain": _validate_domain_param,
    "lat_origin": float,
    "lon_origin": float,
    "datetime_origin": isodate.parse_datetime,
    "forward_duration|backward_duration": isodate.parse_duration,
    "timestep": ("domain_data", isodate.parse_duration),
    "u_vel": _validate_velocity_param,
    "v_vel": _validate_velocity_param,
    "velocity_method": [
        None,
        "lower_troposphere_humidity_weighted",
        "single_height_level",
        "single_pressure_level",
    ],
    "velocity_method_height": [None, float],
    "velocity_method_pressure": [None, float],
}


def build_data_path(root_data_path, trajectory_name):
    data_path = build_data_path_global(
        root_data_path=root_data_path, data_type="trajectory"
    )
    return data_path / "{}.nc".format(trajectory_name)
