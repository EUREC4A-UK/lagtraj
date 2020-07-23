import yaml
import sys
from pathlib import Path

from . import validate_input, build_input_definition_path, examples as input_examples
from .. import DATA_TYPE_PLURAL

FOLDER_STRUCTURE_EXAMPLE = """
data
├── domains
│   ├── eurec4a_circle_eul.yaml
│   └── eurec4a_circle_eul_data
│       ├── an_model_2020-01-01.nc
│       :
│       └── fc_single_2020-01-03.nc
├── forcings
│   ├── eure4a_20191209_12_eul.yaml
│   └── eure4a_20191209_12_eul.nc
└── trajectories
    ├── eure4a_20191209_12_eul.yaml
    └── eure4a_20191209_12_eul.nc
"""


def load_definition(input_name, input_type, root_data_path, required_fields):
    if input_name.startswith("lagtraj://"):
        try:
            input_name = input_name.replace("lagtraj://", "")

            params = input_examples.attempt_read(
                input_name=input_name, input_type=input_type
            )

            input_local_path = build_input_definition_path(
                root_data_path=root_data_path,
                input_name=input_name,
                input_type=input_type,
            )

            if not input_local_path.exists():
                input_local_path.parent.mkdir(exist_ok=True, parents=True)
                with open(input_local_path, "w") as fh:
                    fh.write(yaml.dump(params))

            input_name = input_name
        except input_examples.LagtrajExampleDoesNotExist:
            input_type_plural = DATA_TYPE_PLURAL[input_type]
            print(
                "The requested {} ({}) isn't currently available"
                "in lagtraj\n(maybe you could add it with a pull-request"
                " at https://github.com/EUREC4A-UK/lagtraj/pulls)\n\n"
                "The {} currently defined in lagtraj are:"
                "".format(input_type, input_name, input_type_plural)
            )
            print()
            input_examples.print_available(input_types=[input_type_plural])
            sys.exit(1)
    elif input_name.endswith(".yaml") or input_name.endswith(".yml"):
        # assume we've been passed a full path
        input_local_path = Path(input_name)
        raise Exception(
            "You provided an absolute path to an input"
            " yaml file instead of providing the name"
            " of the input.\nFor for example a trajectory"
            " stored in data/trajectories/eurec4a_20191209_12_lin.yaml"
            " this should be loaded by name as `eurec4a_20191209_12_lin`"
        )
    else:
        input_local_path = build_input_definition_path(
            root_data_path=root_data_path, input_name=input_name, input_type=input_type,
        )
        with open(input_local_path) as fh:
            params = yaml.load(fh, Loader=yaml.FullLoader)

    validate_input(input_params=params, required_fields=required_fields)
    params["name"] = input_name

    return params
