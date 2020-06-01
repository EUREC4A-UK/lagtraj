from pathlib import Path

import xarray as xr
import yaml
import sys

from .. import input_examples
from . import validate_input, build_domain_definition_path


def load_definition(domain, data_path):
    if domain.startswith('lagtraj://'):
        try:
            input_name = domain.replace('lagtraj://', '')

            domain_params = input_examples.attempt_load(
                input_name=input_name, input_type="domain"
            )

            domain_local_path = build_domain_definition_path(
                root_data_path=data_path, domain_name=input_name
            )
            if not domain_local_path.exists():
                domain_local_path.parent.mkdir(exist_ok=True, parents=True)
                with open(domain_local_path, 'w') as fh:
                    fh.write(yaml.dump(domain_params))

            domain_name = input_name
        except input_examples.LagtrajExampleDoesNotExist:
            print("The requested domain ({}) isn't currently available"
                  "in lagtraj\n(maybe you could add it with a pull-request"
                  " at https://github.com/EUREC4A-UK/lagtraj/pulls)\n\n"
                  "The domains currently defined in lagtraj are:"
                  "".format(input_name))
            print()
            input_examples.print_available(input_types=['domains'])
            sys.exit(1)
    else:
        domain_name = domain
        domain_local_path = build_domain_definition_path(
            root_data_path=data_path, domain_name=domain_name
        )
        with open(domain_local_path) as fh:
            domain_params = yaml.load(fh, Loader=yaml.FullLoader)

    validate_input(domain_params=domain_params)
    domain_params['name'] = domain_name

    return domain_params


def load_data(data_path, name):
    domain_data_path = Path(data_path)/name

    ds = xr.open_mfdataset(domain_data_path/"*.nc")

    return ds
