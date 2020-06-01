import os
import datetime
from pathlib import Path

import yaml
import xarray as xr
import numpy as np
import pandas as pd

from .. import DEFAULT_ROOT_DATA_PATH
from .load import load_definition

# Routines for creating a trajectory
# TODO
# - Implement different strategies (single height, weighted over heights, in future possibly hysplit)
# - Add metadata to NetCDF output
# - Improve linear trajectory to work with haversine functions and actual velocities
# - Relax assumption of hourly data?


def create_trajectory():
    pass


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("trajectory_name")
    argparser.add_argument("-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH,
                           type=Path)
    args = argparser.parse_args()

    trajectory_params = load_definition(root_data_path=args.data_path,
                                        trajectory_name=args.trajectory_name)


def get_from_yaml(input_file, directories_file):
    with open(input_file) as this_trajectory_file:
        trajectory_params = yaml.load(this_trajectory_file, Loader=yaml.FullLoader)

    trajectory_type = (trajectory_params["trajectory_type"]).lower()
    if trajectory_type == "eulerian":
        create_eulerian_trajectory(directories_dict, trajectory_params)
    elif trajectory_type == "linear":
        create_linear_trajectory(directories_dict, trajectory_params)
    elif trajectory_type == "single_level":
        create_single_level_trajectory(directories_dict, trajectory_params)
    elif trajectory_type == "weighted":
        create_weighted_trajectory(directories_dict, trajectory_params)
    else:
        raise Exception("Trajectory_type not found")


def create_eulerian_trajectory(directories_dict, trajectory_params):
    times = pd.date_range(
        np.datetime64(trajectory_params["datetime_end"])
        - np.timedelta64(trajectory_params["duration_hours"], "h"),
        periods=trajectory_params["duration_hours"] + 1,
        freq="h",
    )
    nr_hours = len(times)
    lats = np.full((nr_hours), trajectory_params["lat_end"])
    lons = np.full((nr_hours), trajectory_params["lon_end"])
    data = trajectory_to_xarray(times, lats, lons)
    data.to_netcdf(trajectory_filename_parse(directories_dict, trajectory_params))


def create_linear_trajectory(directories_dict, trajectory_params):
    times = pd.date_range(
        np.datetime64(trajectory_params["datetime_end"])
        - np.timedelta64(trajectory_params["duration_hours"], "h"),
        periods=trajectory_params["duration_hours"] + 1,
        freq="h",
    )
    nr_hours = len(times)
    lat_end = trajectory_params["lat_end"]
    lon_end = trajectory_params["lon_end"]
    dlat_dt = trajectory_params["dlat_dt"]
    dlon_dt = trajectory_params["dlon_dt"]
    lats = lat_end - 3600.0 * dlat_dt * np.arange(nr_hours - 1, -1, -1)
    lons = lon_end - 3600.0 * dlon_dt * np.arange(nr_hours - 1, -1, -1)
    data = trajectory_to_xarray(times, lats, lons)
    data.to_netcdf(trajectory_filename_parse(directories_dict, trajectory_params))


def create_single_level_trajectory(directories_dict, trajectory_params):
    raise NotImplementedError


def create_weighted_trajectory(directories_dict, trajectory_params):
    raise NotImplementedError


def trajectory_to_xarray(times, lats, lons):
    lons = (lons) % 360
    data = xr.Dataset({"time": ("time", times)})
    data["lat"] = (("time"), lats)
    data["lon"] = (("time"), lons)
    var_attrs = {
        "lon": {"long_name": "longitude", "units": "degrees_north"},
        "lat": {"long_name": "latitude", "units": "degrees_east"},
    }
    for this_var, this_attr in var_attrs.items():
        data[this_var] = data[this_var].assign_attrs(**this_attr)
    return data


if __name__ == "__main__":
    main()
