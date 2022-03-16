import warnings
from pathlib import Path

import xarray as xr
from tqdm import tqdm

from .. import DEFAULT_ROOT_DATA_PATH
from ..domain.load import load_data as load_domain_data
from ..trajectory.create import (
    required_data_available as required_data_available_for_trajectory,
)
from ..trajectory.load import load_data as load_trajectory_data
from ..utils import optional_debugging, validation
from ..utils.interpolation.levels import make_levels
from ..utils.units import fix_units
from ..utils.xarray import create_attributes_dictionary
from . import build_forcing_data_path, conversion, load, profile_calculation


def _make_latlontime_sampling_points(method, ds_trajectory, ds_domain):
    if method == "domain_data":
        t_min_domn = ds_domain.time.min()
        t_max_domn = ds_domain.time.max()
        t_min_traj = ds_trajectory.time.min()
        t_max_traj = ds_trajectory.time.max()

        if t_min_traj < t_min_domn or t_max_domn < t_max_traj:
            raise Exception(
                "The downloaded domain data does not cover the timespan"
                " of the trajectory requested (ds_trajectory.name)."
                f" t_domain=[{t_min_domn.values},{t_max_domn.values}] and"
                f" t_traj=[{t_min_traj.values},{t_max_traj.values}]. "
                " Please download more domain data"
            )

        da_times = ds_domain.time.sel(time=slice(t_min_traj, t_max_traj))

        ds_sampling = ds_trajectory.interp(
            time=da_times, method="linear", kwargs=dict(bounds_error=True)
        )

        if ds_sampling.dropna(dim="time").time.count() != da_times.count():
            raise Exception(
                "It appears that some of the domain data is incomplete "
                "as interolating the trajectory coordinates to the model "
                "timesteps returned NaNs. Please check that all domain data "
                "has been downloaded."
            )
        # we clip the times to the interval in which the trajectory is defined
        t_min, t_max = ds_trajectory.time.min(), ds_trajectory.time.max()
        ds_sampling = ds_sampling.sel(time=slice(t_min, t_max))
    elif method == "all_trajectory_timesteps":
        ds_sampling = ds_trajectory
    else:
        raise NotImplementedError(
            f"Trajectory sampling method `{method}` not implemented"
        )
    return ds_sampling


def make_forcing(ds_trajectory, ds_domain, levels_definition, sampling_method):
    """
    Make a forcing profiles along ds_trajectory using data in ds_domain.

    See domains.utils.levels.LevelsDefinition and domains.era5.SamplingMethodDefinition
    for how to construct levels_definition and sampling_method objects
    """

    ds_sampling = _make_latlontime_sampling_points(
        method=sampling_method.time_sampling_method,
        ds_trajectory=ds_trajectory,
        ds_domain=ds_domain,
    )

    # `origin_` variables from the trajectory aren't needed for the forcing
    # calculations so lets remove them for now
    ds_sampling = ds_sampling.drop_vars(["origin_lon", "origin_lat", "origin_datetime"])

    ds_sampling["level"] = make_levels(
        method=levels_definition.method,
        n_levels=levels_definition.n_levels,
        z_top=levels_definition.z_top,
        dz_min=levels_definition.dz_min,
    )

    forcing_profiles = []
    for time in tqdm(ds_sampling.time):
        # extract from a single timestep the positions (points in space and
        # time) at which to calculate the forcing profile
        ds_profile_posn = ds_sampling.sel(time=time)
        ds_forcing_profile = profile_calculation.calculate_timestep(
            ds_profile_posn=ds_profile_posn,
            ds_domain=ds_domain,
            sampling_method=sampling_method,
        )
        forcing_profiles.append(ds_forcing_profile)

    ds_forcing = xr.concat(forcing_profiles, dim="time")
    fix_units(ds_forcing)

    ds_forcing["origin_lon"] = ds_trajectory["origin_lon"]
    ds_forcing["origin_lat"] = ds_trajectory["origin_lat"]
    ds_forcing["origin_datetime"] = ds_trajectory["origin_datetime"]

    return ds_forcing


def export(file_path, ds_forcing):
    Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    validation.validate_forcing_profiles(ds_forcing_profiles=ds_forcing)
    encoding = validation.build_valid_encoding(ds=ds_forcing)
    ds_forcing.to_netcdf(file_path, encoding=encoding)


def _validate_existing_forcing(ds_forcing, attr_dict):
    different_attrs = {}
    missing_attrs = {}
    for (k, v) in create_attributes_dictionary(attr_dict).items():
        if k not in ds_forcing.attrs:
            missing_attrs[k] = (None, v)
        elif ds_forcing.attrs[k] != v:
            different_attrs[k] = (ds_forcing.attrs[k], v)

    if len(different_attrs) > 0 or len(missing_attrs) > 0:
        diff_str = "\n".join(
            [
                (
                    f"{k}:\n\tfound: {different_attrs[k][0]}"
                    f"\n\texpected: {different_attrs[k][1]})"
                )
                for k in different_attrs.keys()
            ]
        )
        missing_str = "\n".join(
            [
                (f"{k}: missing\n\texpected: {missing_attrs[k][1]}")
                for k in missing_attrs.keys()
            ]
        )
        raise Exception(f"expected:\n{diff_str}\n{missing_str}")


def main(data_path, forcing_defn, conversion_name=None):
    ds_domain = load_domain_data(root_data_path=data_path, name=forcing_defn.domain)

    try:
        ds_trajectory = load_trajectory_data(
            root_data_path=data_path, name=forcing_defn.trajectory
        )
    except FileNotFoundError:
        raise Exception(
            f"The output file for trajectory `{forcing_defn.trajectory}`"
            " couldn't be found. Please create the trajectory by running: \n"
            f"    python -m lagtraj.trajectory.create {forcing_defn.trajectory}\n"
            "and then run the forcing creation again"
        )

    output_file_path = build_forcing_data_path(
        root_data_path=data_path, forcing_name=forcing_defn.name
    )

    # collect variables which are used for creating attributes of the forcing
    # output file
    attr_dict = dict(
        levels_definition=forcing_defn.levels,
        domain=ds_domain,
        trajectory=ds_trajectory,
        sampling_method=forcing_defn.sampling,
        trajectory_name=forcing_defn.trajectory,
    )

    if output_file_path.exists():
        ds_forcing = xr.open_dataset(output_file_path)
        if conversion_name is None:
            raise Exception(
                f"A file already exists at the path `{output_file_path}`. "
                "Please delete this file and run the forcing creation again"
            )
        else:
            try:
                _validate_existing_forcing(ds_forcing=ds_forcing, attr_dict=attr_dict)
                warnings.warn(
                    f"Using existing forcing file `{output_file_path}` to convert to "
                    f"`{conversion_name}` format."
                )
            except Exception as e:
                raise Exception(
                    f"A forcing file already exists at the path `{output_file_path}`. "
                    f"Applying the `{conversion_name}` conversion was halted "
                    "because the following attributes of the file were different than "
                    f"expected:\n{e}"
                )
    else:
        ds_forcing = make_forcing(
            levels_definition=forcing_defn.levels,
            ds_domain=ds_domain,
            ds_trajectory=ds_trajectory,
            sampling_method=forcing_defn.sampling,
        )
        ds_forcing.attrs.update(create_attributes_dictionary(attr_dict))

    if not output_file_path.exists():
        export(
            ds_forcing=ds_forcing,
            file_path=output_file_path,
        )
        print("Wrote forcing file to `{}`".format(output_file_path))

    if conversion_name is not None:
        converted_output_file_path = conversion.export_for_target(
            ds_forcing=ds_forcing,
            conversion_name=conversion_name,
            root_data_path=data_path,
        )
        print("Wrote converted forcing file to `{}`".format(converted_output_file_path))


def _make_cli_argparser():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("forcing")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    available_conversion_targets = conversion.targets.available.keys()
    argparser.add_argument(
        "-c",
        "--conversion",
        help="name of output conversion to use, available conversions: "
        f"{', '.join(available_conversion_targets)}",
        default=None,
    )
    argparser.add_argument("--debug", default=False, action="store_true")
    return argparser


def cli(args=None):
    """
    Function called with arguments passed from the command line when making
    trajectories through the CLI. When `args==None` they will be taken from
    `sys.argv`
    """
    argparser = _make_cli_argparser()
    args = argparser.parse_args(args=args)

    forcing_defn = load.load_definition(
        root_data_path=args.data_path, forcing_name=args.forcing
    )
    with optional_debugging(args.debug):
        main(
            data_path=args.data_path,
            forcing_defn=forcing_defn,
            conversion_name=args.conversion,
        )


def has_data_for_cli_command(args):
    argparser = _make_cli_argparser()
    args = argparser.parse_args(args=args)

    data_path = args.data_path
    forcing_defn = load.load_definition(
        root_data_path=args.data_path, forcing_name=args.forcing
    )
    trajectory_name = forcing_defn.trajectory

    return required_data_available_for_trajectory(
        data_path=data_path, trajectory_name=trajectory_name
    )


if __name__ == "__main__":
    cli()
