import semver


from .. import build_data_path


def build_input_definition_path(root_data_path, input_name, input_type):
    data_path = build_data_path(root_data_path=root_data_path, data_type=input_type)

    return data_path / (input_name + ".yaml")


class InvalidInputDefinition(Exception):
    pass


def validate_input(input_params, required_fields):
    def _check_field(f_name, f_option):
        # allows for an optional parameter by putting `None` in the list of
        # options
        missing_allowed = type(f_option) in [list, tuple] and None in f_option

        if missing_allowed and f_name not in input_params:
            pass
        elif f_name not in input_params:
            raise InvalidInputDefinition("Missing `{}` field".format(f_name))
        elif type(f_option) == type and type(input_params[f_name]) != f_option:
            raise InvalidInputDefinition(
                "Field `{}` should have type {}, but has type `{}`".format(
                    f_name, f_option, type(input_params[f_name])
                )
            )
        elif callable(f_option):
            return f_option(input_params[f_name])
        elif (
            type(f_option) == str
            or type(f_option) == int
            or type(f_option) == float
            and input_params[f_name] == f_option
        ):
            # allows setting an option to be just a string
            pass
        elif type(f_option) == list or type(f_option) == tuple:
            f_options = f_option
            exceptions = []
            new_val = None
            for f_option in f_options:
                try:
                    new_val = _check_field(f_name, f_option)
                except Exception as e:
                    exceptions.append(e)
            if len(exceptions) == len(f_options):
                raise InvalidInputDefinition(exceptions)
            else:
                if new_val is not None:
                    return new_val
        else:
            raise NotImplementedError(f_name, f_option)

    checked_valid_fields = []
    for f_name, f_option in required_fields.items():
        if "|" in f_name:
            f_names = f_name.split("|")
            exceptions = []
            for f_name in f_names:
                try:
                    new_val = _check_field(f_name, f_option)
                    if new_val is not None:
                        input_params[f_name] = new_val
                except Exception as e:
                    exceptions.append(e)

            if len(exceptions) == len(f_names):
                raise InvalidInputDefinition(
                    "None of the fields `{}` were correctly defined in the input"
                    " definition. Please define at least one with the {} type"
                    "\n\nThe errors incurred where: {}"
                    "".format(", ".join(f_names), f_option, exceptions)
                )
            else:
                checked_valid_fields += f_names
        else:
            new_val = _check_field(f_name, f_option)
            if new_val is not None:
                input_params[f_name] = new_val
            checked_valid_fields.append(f_name)

    if "version" in input_params:
        try:
            semver.parse(str(input_params["version"]))
        except ValueError:
            raise InvalidInputDefinition(
                "Versioning labels should follow the semver convention "
                "(http://semver.org) of `MINOR.MAJOR.PATCH` (e.g. the "
                "first version you make might be `1.0.0`)."
                f" The version is currenty given as `{input_params['version']}`."
            )
        else:
            checked_valid_fields.append("version")

    extra_fields = set(input_params.keys()).difference(checked_valid_fields)

    if len(extra_fields) > 0:
        raise InvalidInputDefinition(
            "Input definition has the following"
            " extra fields: {}".format(", ".join(extra_fields))
        )
