from pathlib import Path

import xarray as xr
from tqdm import tqdm

from .. import DEFAULT_ROOT_DATA_PATH
from . import era5, load, build_forcing_data_path
from .utils.levels import make_levels
from ..utils import optional_debugging
from ..domain.load import load_data as load_domain_data
from ..trajectory.load import load_data as load_trajectory_data


def _make_latlontime_sampling_points(method, ds_trajectory, ds_domain):
    if method == "model_timesteps":
        da_time = ds_domain.time
        da_sampling = ds_trajectory.interp(time=da_time, method="linear")
        # we clip the times to the interval in which the trajectory is defined
        t_min, t_max = ds_trajectory.time.min(), ds_trajectory.time.max()
        da_sampling = da_sampling.sel(time=slice(t_min, t_max))
    elif method == "all_trajectory_timesteps":
        da_sampling = ds_trajectory
    else:
        raise NotImplementedError(
            "Trajectory sampling method `{}` not implemented".format()
        )
    return da_sampling


def make_forcing(
    ds_trajectory, ds_domain, levels_definition, sampling_method,
):
    """
    Make a forcing profiles along ds_trajectory using data in ds_domain.

    See domains.utils.levels.LevelsDefinition and domains.era5.SamplingMethodDefinition
    for how to construct levels_definition and sampling_method objects
    """

    da_sampling = _make_latlontime_sampling_points(
        method=sampling_method.time_sampling_method,
        ds_trajectory=ds_trajectory,
        ds_domain=ds_domain,
    )

    da_sampling["levels"] = make_levels(
        method=levels_definition.method,
        n_levels=levels_definition.n_levels,
        z_top=levels_definition.z_top,
        dz_min=levels_definition.dz_min,
    )

    ds_forcing = xr.Dataset(coords=da_sampling.coords,)

    # XXX: eventually this will be replaced by a function which doesn't assume
    # that the domain data is era5 data
    timestep_function = era5.calculate_timestep

    ds_timesteps = []
    for time in tqdm(ds_forcing.time):
        da_pt = da_sampling.sel(time=time)
        ds_timestep = timestep_function(
            da_pt=da_pt, ds_domain=ds_domain, sampling_method=sampling_method
        )
        ds_timesteps.append(ds_timestep)

    return xr.concat(ds_timesteps, dim="time")


def export(file_path, ds_forcing, format):
    if format != "dummy_netcdf":
        raise NotImplementedError(format)
    # TODO: add hightune format export here

    Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    ds_forcing.to_netcdf(file_path)


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("forcing")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument("-f", "--output-format", default="dummy_netcdf", type=str)
    argparser.add_argument("--debug", default=False, action="store_true")
    args = argparser.parse_args()

    forcing_defn = load.load_definition(
        root_data_path=args.data_path, forcing_name=args.forcing
    )

    ds_domain = load_domain_data(
        root_data_path=args.data_path, name=forcing_defn.domain
    )
    ds_trajectory = load_trajectory_data(
        root_data_path=args.data_path, name=forcing_defn.trajectory
    )

    with optional_debugging(args.debug):
        ds_forcing = make_forcing(
            levels_definition=forcing_defn.levels,
            ds_domain=ds_domain,
            ds_trajectory=ds_trajectory,
            sampling_method=forcing_defn.sampling,
        )

    output_file_path = build_forcing_data_path(
        root_data_path=args.data_path, forcing_name=forcing_defn.name
    )

    export(ds_forcing=ds_forcing, file_path=output_file_path, format=args.output_format)

    print("Wrote forcing file to `{}`".format(output_file_path))


if __name__ == "__main__":
    main()
