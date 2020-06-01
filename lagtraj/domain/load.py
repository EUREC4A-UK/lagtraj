from ..input_definitions import load
from . import INPUT_REQUIRED_FIELDS, build_domain_data_path
from .sources import era5


def load_definition(root_data_path, name):
    return load.load_definition(
        root_data_path=root_data_path,
        input_name=name,
        input_type="domain",
        required_fields=INPUT_REQUIRED_FIELDS,
    )


def load_data(root_data_path, name):
    # check that the domain data definition exists first, might be needed to
    # read the data in
    try:
        domain_def = load_definition(root_data_path=root_data_path, name=name)
    except FileNotFoundError:
        print("The domain definition for `{}` couldn't be found".format(name))
        raise

    data_path = build_domain_data_path(
        root_data_path=root_data_path, domain_name=domain_def["name"]
    )
    if domain_def["source"] == "era5":
        ds = era5.load_data(data_path=data_path)
    else:
        raise NotImplementedError(domain_def["source"])

    return ds
