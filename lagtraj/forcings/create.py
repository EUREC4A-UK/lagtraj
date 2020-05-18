import yaml
import datetime
import xarray as xr
import os
import numpy as np
import pandas as pd
from lagtraj.utils.parsers import (
    domain_filename_parse,
    trajectory_filename_parse,
    forcings_filename_parse,
)
from lagtraj.utils.levels import make_levels

# from lagtraj.utils.era5 import add_heights_and_pressures
from lagtraj.utils.hightune import hightune_variables


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
    with open(input_file) as this_forcings_file:
        forcings_dict = yaml.load(this_forcings_file, Loader=yaml.FullLoader)

    trajectory_file = trajectory_filename_parse(
        directories_dict, forcings_dict["trajectory"]
    )
    trajectory_data = xr.open_dataset(trajectory_file)
    times = trajectory_data["time"]

    levels = make_levels(forcings_dict)

    forcing_data = xr.Dataset()
    for timestep in trajectory_data["time"]:
        append_timestep(forcings_dict, forcing_data, trajectory_data, timestep)
    export_to_hightune(forcing_data)


def append_timestep(forcings_dict, forcing_data, trajectory_data, timestep):
    if (forcings_dict["source"]).lower() == "era5":
        append_era5_timestep(forcings_dict, forcing_data, trajectory_data, timestep)


def append_era5_timestep(forcings_dict, forcing_data, trajectory_data, timestep):
    # FOR EACH TIME STEP
    # find the right input data
    # reinterpolate data to height or pressure with effective height (tbd)
    # create 'mask' for forcings on domain
    # calculate profiles and forcings
    lon_indices, lat_indices = calculate_current_lat_lon_range(
        forcings_dict, trajectory_data, timestep
    )
    timestep_input_data = xr.Dataset()
    append_single_level_an(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    append_single_level_fc(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    append_model_level_an(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    append_model_level_fc(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    add_heights_and_pressures(timestep_input_data)
    timestep_interpolated_data = interpolate_to_height_levels(timestep_input_data)
    timestep_forcing_data = calculate_forcings(
        forcings_dict, timestep_interpolated_data
    )
    # append to forcing_data using xarray concatenation


def calculate_current_lat_lon_range(forcings_dict, trajectory_data, timestep):
    lon_indices = []
    lat_indices = []
    return lon_indices, lat_indices


def append_single_level_an(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def append_single_level_fc(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def append_model_level_an(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def append_model_level_fc(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def add_heights_and_pressures(timestep_input_data):
    pass


def interpolate_to_height_levels(timestep_input_data):
    timestep_interpolated_data = []
    return timestep_interpolated_data


def calculate_forcings(forcings_dict, forcing_interpolated_data):
    timestep_forcing_data = []
    return timestep_forcing_data


def export_to_hightune(forcing_data):
    pass


def export_to_ecmwf(forcing_data):
    pass


if __name__ == "__main__":
    main()
