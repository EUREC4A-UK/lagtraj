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

from . import era5

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
    """
    forcing.yaml:
        trajectory_name: ...
        time_sampling_method: "model_timesteps"
        domain: eul_...
        source: era5
    """

    with open(directories_file) as this_directories_file:
        directories_dict = yaml.load(this_directories_file, Loader=yaml.FullLoader)
    with open(input_file) as this_forcings_file:
        forcings_dict = yaml.load(this_forcings_file, Loader=yaml.FullLoader)

    trajectory_file = trajectory_filename_parse(
        directories_dict, forcings_dict["trajectory"]
    )
    ds_trajectory = xr.open_dataset(trajectory_file)

    time_sampling_method = trajectory_file.get(
        'time_sampling_method', default="model_timesteps"
    )

    if time_sampling_method == "model_timesteps":
        if forcings_dict['source'] == "era5":
            da_time = era5.get_available_timesteps(
                domain=forcings_dict['source']
            )
        else:
            raise NotImplementedError(forcings_dict['source'])

        da_sampling = ds_trajectory.resample(time=da_time).interpolate("linear")
    elif time_sampling_method == "all_trajectory_timesteps":
        da_sampling = ds_trajectory
    else:
        raise NotImplementedError(
            "Trajectory sampling method `{}` not implemented".format(
            )
        )

    da_sampling['levels'] = make_levels(forcings_dict)

    ds_forcings = xr.Dataset(
        coords=da_sampling.coords,
    )

    if forcings_dict['source'] == "era5":
        timestep_function = era5.calculate_timestep
    else:
        raise NotImplementedError(forcings_dict['source'])

    ds_timesteps = []
    for time in forcing_data.time:
        ds_timestep = timestep_function(time, domain=forcings_dict['domain'])
        ds_timesteps.append(ds_timestep)

    ds_forcings = xr.concat(ds_timesteps, dim='time')

    export_to_hightune(ds_forcings)
