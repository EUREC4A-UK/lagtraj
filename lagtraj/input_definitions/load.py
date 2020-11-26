import yaml
import sys
from pathlib import Path

from . import validate_input, build_input_definition_path, examples as input_examples
from .. import DATA_TYPE_PLURAL
from .examples import get_available as get_available_input_examples

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
                "The requested {} ({}) isn't currently available "
                "in lagtraj\n(maybe you could add it with a pull-request"
                " at https://github.com/EUREC4A-UK/lagtraj/pulls)\n\n"
                "The {} currently defined in lagtraj are:"
                "".format(input_type, input_name, input_type_plural)
            )
            print()
            input_examples.print_available(input_types=[input_type_plural])
            sys.exit(1)
    else:
        if input_name.endswith(".yaml") or input_name.endswith(".yml"):
            # assume we've been passed a full path
            input_local_path = Path(input_name)
            if not input_local_path.exists():
                raise Exception(
                    "You provided an absolute path to an input"
                    " yaml file, but that file doesn't appear"
                    " to exist. Maybe you intended to just"
                    " load this input definition by name instead?"
                )
            # check that this provided yaml-file is in the correct folder
            # structure
            input_type_plural = DATA_TYPE_PLURAL[input_type]
            if not (
                input_local_path.parent.parent.name == "data"
                and input_local_path.parent.name == input_type_plural
            ):
                lagtraj_input_examples = list(
                    get_available_input_examples(input_types=[input_type])
                )
                s = ", ".join(lagtraj_input_examples[:3])  #  show the first three only
                raise Exception(
                    "The yaml input-file you provided does not"
                    " exist in the correct direction structure."
                    " lagtraj assumes that all data is stored in"
                    " a directory structure as follows (so that"
                    " so that the relevant input data and output"
                    f" can be found by lagtraj): \n{FOLDER_STRUCTURE_EXAMPLE}"
                    f"{input_type_plural} input definitions bundled with `lagtraj` can be used as well (for example {s}"
                    " - see all available with `python -m lagtraj.input_definitions.examples`)"
                    "e.g. \n"
                    "python -m lagtraj.trajectory.create lagtraj://eurec4a_20200202_12_lag \n"
                    "The input files are then copied over to the data directory."
                )
        else:
            input_local_path = build_input_definition_path(
                root_data_path=root_data_path,
                input_name=input_name,
                input_type=input_type,
            )
            if not input_local_path.exists():
                input_type_plural = DATA_TYPE_PLURAL[input_type]
                raise Exception(
                    f"The requested {input_type} ({input_name}) wasn't found. "
                    f"To use a {input_type} with this name please define its "
                    f"parameters {input_local_path}\n"
                    "(or run `python -m lagtraj.input_definitions.examples` all available with"
                    "to see the ones currently bundled with lagtraj"
                )
            print()
        with open(input_local_path) as fh:
            params = yaml.load(fh, Loader=yaml.FullLoader)

    validate_input(input_params=params, required_fields=required_fields)
    params["name"] = input_name

    if not "version" in params:
        params["version"] = "unversioned"

    return params
