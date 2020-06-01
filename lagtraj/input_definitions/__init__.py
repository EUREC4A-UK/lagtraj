from .. import build_data_path


def build_input_definition_path(root_data_path, input_name, input_type):
    data_path = build_data_path(root_data_path=root_data_path,
                                data_type=input_type)

    return data_path / (input_name + ".yaml")


class InvalidInputDefinition(Exception):
    pass


def validate_input(input_params, required_fields):
    def _check_field(f_name, f_type):
        if f_name not in input_params:
            raise InvalidInputDefinition("Missing `{}` field".format(f_name))
        if callable(f_type):
            try:
                f_type(input_params[f_name])
            except Exception:
                print("There was an issue parsing the `{}` field from the input"
                      " definition".format(f_name))
                raise
        elif type(input_params[f_name]) != f_type:
            raise InvalidInputDefinition(
                "Field `{}` should have type {}, but has type `{}`".format(
                 f_name, f_type, type(input_params[f_name]))
            )

    checked_valid_fields = []
    for f_name, f_type in required_fields.items():
        if "|" in f_name:
            f_names = f_name.split("|")
            exceptions = []
            for f_name in f_names:
                try:
                    _check_field(f_name, f_type)
                except Exception as e:
                    exceptions.append(e)

            if len(exceptions) == len(f_names):
                raise InvalidInputDefinition(
                    "None of the fields `{}` were correctly defined in the input"
                    " definition. Please define at least one with the {} type"
                    "".format(", ".join(f_names), f_type))
            else:
                checked_valid_fields += f_names
        else:
            _check_field(f_name, f_type)
            checked_valid_fields.append(f_name)

    extra_fields = set(input_params.keys()).difference(checked_valid_fields)

    if len(extra_fields) > 0:
        raise InvalidInputDefinition(
            "Input definition has the following"
            " extra fields: {}".format(", ".join(extra_fields))
        )
