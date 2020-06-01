from collections import namedtuple
from pathlib import Path


LatLonBoundingBox = namedtuple("LatLonBoundingBox",
                               ["lat_min", "lat_max", "lon_min", "lon_max"]
                               )

LatLonSamplingResolution = namedtuple("LatLonSamplingResolution",
                                      ["lat", "lon"]
                                      )


INPUT_REQUIRED_FIELDS = dict(
    source=str,
    lat_min=float,
    lat_max=float,
    lon_min=float,
    lon_max=float,
    lat_samp=float,
    lon_samp=float
)


class InvalidDomainDefinition(Exception):
    pass


def validate_input(domain_params):
    for f, f_type in INPUT_REQUIRED_FIELDS.items():
        if f not in domain_params:
            raise InvalidDomainDefinition("Missing `{}` field".format(f))
        elif type(domain_params[f]) != f_type:
            raise InvalidDomainDefinition("Field `{}` should have type {}"
                                          "".format(f, f_type))

    extra_fields = set(
        domain_params.keys()
    ).difference(INPUT_REQUIRED_FIELDS.keys())

    if len(extra_fields) > 0:
        raise InvalidDomainDefinition("Domain definition has the following"
                                      " extra fields: {}".format(
                                          ", ".join(extra_fields)
                                      ))


def build_domain_data_path(root_data_path, domain_name):
    return Path(root_data_path)/"domains"/(domain_name + "_data")


def build_domain_definition_path(root_data_path, domain_name):
    return Path(root_data_path)/"domains"/(domain_name + ".yaml")
