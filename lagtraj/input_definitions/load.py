import yaml
import sys

from . import validate_input, build_input_definition_path, examples as input_examples
from .. import DATA_TYPE_PLURAL


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
    else:
        input_local_path = build_input_definition_path(
            root_data_path=root_data_path, input_name=input_name, input_type=input_type,
        )
        with open(input_local_path) as fh:
            params = yaml.load(fh, Loader=yaml.FullLoader)

    validate_input(input_params=params, required_fields=required_fields)
    params["name"] = input_name

    return params
