from pathlib import Path

import xarray as xr

from ..input_definitions import load
from . import INPUT_REQUIRED_FIELDS
from .. import build_data_path


def load_definition(data_path, domain_name):
    return load.load_definition(data_path=data_path, input_name=domain_name,
                                input_type="domain",
                                required_fields=INPUT_REQUIRED_FIELDS)


def load_data(root_data_path, name):
    data_path = build_data_path(root_data_path=root_data_path,
                                data_type="domain")

    domain_data_path = data_path/"{}_data/*.nc".format(name)
    ds = xr.open_mfdataset(domain_data_path)

    return ds
