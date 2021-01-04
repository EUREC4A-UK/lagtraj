import warnings
from pathlib import Path
from tqdm import tqdm
import xarray as xr

from .. import DEFAULT_ROOT_DATA_PATH
from . import profile_calculation, load, build_forcing_data_path, conversion
from ..utils.interpolation.levels import make_levels
from ..utils import optional_debugging, validation
from ..domain.load import load_data as load_domain_data
from ..trajectory.load import load_data as load_trajectory_data
from ..utils.xarray import create_attributes_dictionary
from ..utils.units import fix_units


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
    return ds_forcing


def export(file_path, ds_forcing, format):
    if format != "dummy_netcdf":
        raise NotImplementedError(format)
    # TODO: add dephy format export here

    Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    validation.validate_forcing_profiles(ds_forcing_profiles=ds_forcing)
    encoding = validation.build_valid_encoding(ds=ds_forcing)
    ds_forcing.to_netcdf(file_path, encoding=encoding)


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("forcing")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument("-f", "--output-format", default="dummy_netcdf", type=str)
    available_conversion_targets = conversion.targets.available.keys()
    argparser.add_argument(
        "-t",
        "--target-model",
        help="name of the model to target, available targets: "
        f"{', '.join(available_conversion_targets)}",
        default=None,
    )
    argparser.add_argument("--debug", default=False, action="store_true")
    args = argparser.parse_args()

    forcing_defn = load.load_definition(
        root_data_path=args.data_path, forcing_name=args.forcing
    )

    with optional_debugging(args.debug):
        ds_domain = load_domain_data(
            root_data_path=args.data_path, name=forcing_defn.domain
        )
    try:
        ds_trajectory = load_trajectory_data(
            root_data_path=args.data_path, name=forcing_defn.trajectory
        )
    except FileNotFoundError:
        raise Exception(
            f"The output file for trajectory `{forcing_defn.trajectory}`"
            " couldn't be found. Please create the trajectory by running: \n"
            f"    python -m lagtraj.trajectory.create {forcing_defn.trajectory}\n"
            "and then run the forcing creation again"
        )

    output_file_path = build_forcing_data_path(
        root_data_path=args.data_path, forcing_name=forcing_defn.name
    )

    # collect variables which are used for creating attributes of the forcing
    # output file
    attr_dict = dict(
        levels_definition=forcing_defn.levels,
        ds_domain=ds_domain,
        ds_trajectory=ds_trajectory,
        sampling_method=forcing_defn.sampling,
        trajectory_name=forcing_defn.trajectory,
    )

    if output_file_path.exists():
        ds_forcing = xr.open_dataset(output_file_path)
        if args.target_model is None:
            raise Exception(
                f"A file already exists at the path `{output_file_path}`. "
                "Please delete this file and run the forcing creation again"
            )
        else:
            with optional_debugging(args.debug):
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
                            f"{k}:\n\tfound: {different_attrs[k][0]}\n\texpected: {different_attrs[k][1]})"
                        )
                        for k in different_attrs.keys()
                    ]
                )
                missing_str = "\n".join(
                    [
                        (f"{k}: missing\n\texpected: {missing_attrs[k][1]})")
                        for k in missing_attrs.keys()
                    ]
                )
                raise Exception(
                    f"A forcing file already exists at the path `{output_file_path}`. "
                    f"Applying the `{args.target_model}` conversion was halted "
                    "because the following attributes of the file were different than "
                    f"expected:\n{diff_str}\n{missing_str}"
                )
            else:
                warnings.warn(
                    f"Using existing forcing file `{output_file_path}` to convert to "
                    f"`{args.target_model}` format."
                )

    else:
        with optional_debugging(args.debug):
            ds_forcing = make_forcing(
                levels_definition=forcing_defn.levels,
                ds_domain=ds_domain,
                ds_trajectory=ds_trajectory,
                sampling_method=forcing_defn.sampling,
            )

        ds_forcing["origin_lon"] = ds_trajectory["origin_lon"]
        ds_forcing["origin_lat"] = ds_trajectory["origin_lat"]
        ds_forcing["origin_datetime"] = ds_trajectory["origin_datetime"]
        ds_forcing.attrs.update(ds_domain.attrs)
        ds_forcing.attrs.update(create_attributes_dictionary(attr_dict))

    if args.target_model is not None:
        conversion.export_for_target(
            ds_forcing=ds_forcing,
            target_name=args.target_model,
            root_data_path=args.data_path,
        )
    else:
        export(
            ds_forcing=ds_forcing, file_path=output_file_path, format=args.output_format
        )
        print("Wrote forcing file to `{}`".format(output_file_path))


if __name__ == "__main__":
    main()
