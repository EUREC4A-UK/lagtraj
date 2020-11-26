from ..input_definitions import load
from . import INPUT_REQUIRED_FIELDS, build_domain_data_path
from .sources import era5


def load_definition(domain_name, data_path):
    domain_params = load.load_definition(
        input_name=domain_name,
        input_type="domain",
        root_data_path=data_path,
        required_fields=INPUT_REQUIRED_FIELDS,
    )

    if not "version" in domain_params:
        domain_params["version"] = "unversioned"

    return domain_params


def load_data(root_data_path, name):
    # check that the domain data definition exists first, might be needed to
    # read the data in
    try:
        domain_def = load_definition(data_path=root_data_path, domain_name=name)
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

    ds.attrs["data_source"] = domain_def["source"]

    if ds.version != domain_def["version"]:
        raise Exception("The domain data in `{data_path}` doesn't match the version"
                        " stored in the data definition for `{name}`. Please delete"
                        " the domain data and re-download.")

    return ds
