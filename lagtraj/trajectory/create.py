import yaml
import datetime
import xarray as xr
import os
import numpy as np
import pandas as pd
from lagtraj.utils.parsers import domain_filename_parse, trajectory_filename_parse

""" Routines for creating a trajectory

Todo
- Implement different strategies (single height, weighted over heights, in future possibly hysplit)
- Add metadata to NetCDF output
- Improve linear trajectory to work with haversine functions and actual velocities
- Relax assumption of hourly data?"""


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("input_file")
    argparser.add_argument(
        "-f", "--file", dest="directories_file", default="directories.yaml", type=str
    )
    args = argparser.parse_args()
    get_from_yaml(args.input_file, args.directories_file)


def get_from_yaml(input_file, directories_file):
    with open(directories_file) as this_directories_file:
        directories_dict = yaml.load(this_directories_file, Loader=yaml.FullLoader)
    with open(input_file) as this_trajectory_file:
        trajectory_dict = yaml.load(this_trajectory_file, Loader=yaml.FullLoader)

    trajectory_type = (trajectory_dict["trajectory_type"]).lower()
    if trajectory_type == "eulerian":
        create_eulerian_trajectory(directories_dict, trajectory_dict)
    elif trajectory_type == "linear":
        create_linear_trajectory(directories_dict, trajectory_dict)
    elif trajectory_type == "single_level":
        create_single_level_trajectory(directories_dict, trajectory_dict)
    elif trajectory_type == "weighted":
        create_weighted_trajectory(directories_dict, trajectory_dict)
    else:
        raise Exception("Trajectory_type not found")


def create_eulerian_trajectory(directories_dict, trajectory_dict):
    times = pd.date_range(
        np.datetime64(trajectory_dict["datetime_end"])
        - np.timedelta64(trajectory_dict["duration_hours"], "h"),
        periods=trajectory_dict["duration_hours"] + 1,
        freq="h",
    )
    nr_hours = len(times)
    lats = np.full((nr_hours), trajectory_dict["lat_end"])
    lons = np.full((nr_hours), trajectory_dict["lon_end"])
    data = trajectory_to_xarray(times, lats, lons)
    data.to_netcdf(trajectory_filename_parse(directories_dict, trajectory_dict))


def create_linear_trajectory(directories_dict, trajectory_dict):
    times = pd.date_range(
        np.datetime64(trajectory_dict["datetime_end"])
        - np.timedelta64(trajectory_dict["duration_hours"], "h"),
        periods=trajectory_dict["duration_hours"] + 1,
        freq="h",
    )
    nr_hours = len(times)
    lat_end = trajectory_dict["lat_end"]
    lon_end = trajectory_dict["lon_end"]
    dlat_dt = trajectory_dict["dlat_dt"]
    dlon_dt = trajectory_dict["dlon_dt"]
    lats = lat_end - 3600.0 * dlat_dt * np.arange(nr_hours - 1, -1, -1)
    lons = lon_end - 3600.0 * dlon_dt * np.arange(nr_hours - 1, -1, -1)
    data = trajectory_to_xarray(times, lats, lons)
    data.to_netcdf(trajectory_filename_parse(directories_dict, trajectory_dict))


def create_single_level_trajectory(directories_dict, trajectory_dict):
    pass


def create_weighted_trajectory(directories_dict, trajectory_dict):
    pass


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
