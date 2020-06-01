import xarray as xr

from ..input_definitions import load
from . import INPUT_REQUIRED_FIELDS, build_data_path


def load_definition(root_data_path, trajectory_name):
    return load.load_definition(root_data_path=root_data_path,
                                input_name=trajectory_name,
                                input_type="trajectory",
                                required_fields=INPUT_REQUIRED_FIELDS)


def load_data(root_data_path, name):
    data_path = build_data_path(root_data_path=root_data_path,
                                data_type="trajectory")
    trajectory_data_path = data_path/"{}.nc".format(name)

    ds = xr.open_dataset(trajectory_data_path)
    return ds
