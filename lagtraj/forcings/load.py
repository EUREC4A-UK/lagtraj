import xarray as xr

from ..input_definitions import load
from . import (
    INPUT_REQUIRED_FIELDS,
    ForcingDefinition,
    ForcingLevelsDefinition,
    ForcingSamplingDefinition,
    build_forcing_data_path,
)


def load_definition(root_data_path, forcing_name):
    forcing_params = load.load_definition(
        root_data_path=root_data_path,
        input_name=forcing_name,
        input_type="forcing",
        required_fields=INPUT_REQUIRED_FIELDS,
    )

    forcing_levels_definition = ForcingLevelsDefinition(
        n_levels=forcing_params["levels_number"],
        z_top=forcing_params["levels_ztop"],
        dz_min=forcing_params.get("levels_dzmin", None),
        method=forcing_params["levels_method"],
    )

    sampling_definition = ForcingSamplingDefinition(
        gradient_method=forcing_params["gradient_method"],
        advection_velocity_sampling_method=forcing_params[
            "advection_velocity_sampling_method"
        ],
        averaging_width=forcing_params["averaging_width"],
        time_sampling_method=forcing_params.get("time_sampling_method", "domain_data"),
        mask=forcing_params.get("sampling_mask"),
    )

    forcing_definition = ForcingDefinition(
        trajectory=forcing_params["trajectory"],
        domain=forcing_params["domain"],
        sampling=sampling_definition,
        levels=forcing_levels_definition,
        name=forcing_params["name"],
        version=forcing_params["version"],
    )

    return forcing_definition


def load_data(root_data_path, forcing_name):
    """Function to load the data after creation"""
    forcing_definition = load_definition(
        root_data_path=root_data_path, forcing_name=forcing_name
    )

    forcing_data_path = build_forcing_data_path(
        root_data_path=root_data_path, forcing_name=forcing_name
    )
    ds = xr.open_dataset(forcing_data_path)

    if "version" not in ds.attrs:
        raise Exception(
            f"The forcing stored in `{forcing_data_path}` "
            "doesn't `version` attribute set. Please delete and "
            "recreate the forcing or set the correct version. "
            "Was expecting to find forcing with version "
            f"`{forcing_definition['version']}`"
        )

    if ds.version != forcing_definition.version:
        raise Exception(
            f"The version of the forcing stored in `{forcing_data_path}` "
            "doesn't match the version in the input definition yaml file "
            f"for a forcing named `{forcing_name}`. Please delete the forcing "
            "netCDF file and recreate it"
        )

    return ds
