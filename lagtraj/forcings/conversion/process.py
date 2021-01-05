"""Module that enables conversions from era5 format to other formats
"""

from pathlib import Path

from ...utils import validation
from ...utils.interpolation.levels import make_levels
from ... import DEFAULT_ROOT_DATA_PATH
from .targets import available as available_targets
from .load import load_definition as load_conversion_defn
from .. import build_forcing_data_path


def export(ds_forcing, output_filepath, conversion_defn):
    if conversion_defn.levels.method is None or conversion_defn.levels.method == "copy":
        da_levels = ds_forcing["level"]
    else:
        da_levels = make_levels(
            method=conversion_defn.levels.method,
            n_levels=conversion_defn.levels.n_levels,
            z_top=conversion_defn.levels.z_top,
            dz_min=conversion_defn.levels.dz_min,
        )

    if conversion_defn.export_format not in available_targets:
        raise Exception(
            "The export format `{conversion_defn.export_format}` "
            "isn't currently supported. Target models currently "
            "supported: {', '.join(available_targets.keys())}."
        )

    conversion_module = available_targets[conversion_defn.export_format]
    conversion_func = getattr(conversion_module, "from_era5")
    if conversion_func is None:
        raise Exception(
            "The module describing conversion to the "
            "`{conversion_defn.export_format}` ({conversion_module.__name__})"
        )

    ds_converted = conversion_func(
        ds_forcing,
        da_levels,
        parameters=conversion_defn.parameters,
        metadata=conversion_defn.metadata,
    )

    Path(output_filepath).parent.mkdir(parents=True, exist_ok=True)

    # set up the correct encoding for the output (ensures correct units on time
    # coordinates etc) and check if the target model has some specific
    # requirements for how the netCDF file is written too
    export_kwargs = dict(encoding=validation.build_valid_encoding(ds=ds_converted))
    export_kwargs.update(conversion_module.getattr("EXPORT_KWARGS", {}))

    ds_converted.to_netcdf(output_filepath, **export_kwargs)


def export_for_target(ds_forcing, target_name, root_data_path=DEFAULT_ROOT_DATA_PATH):
    """
    Export the forcing `ds_forcing` into `root_data_path` by applying
    conversion identified by `target_name`
    """
    if "name" not in ds_forcing.attrs:
        raise Exception(
            "To be able to export a forcing to file its `name` attribute"
            " must be set."
        )

    conversion_defn = load_conversion_defn(
        root_data_path=root_data_path,
        forcing_name=ds_forcing.name,
        target_name=target_name,
    )
    output_filepath = build_forcing_data_path(
        root_data_path=root_data_path,
        forcing_name=ds_forcing.name,
        target_name=target_name,
    )
    export(
        ds_forcing=ds_forcing,
        output_filepath=output_filepath,
        conversion_defn=conversion_defn,
    )
