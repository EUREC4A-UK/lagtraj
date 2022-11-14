import semver

from .. import build_data_path


def build_input_definition_path(
    root_data_path, input_name, input_type, input_subtype=None
):
    data_path = build_data_path(root_data_path=root_data_path, data_type=input_type)

    if input_subtype is None:
        filename = f"{input_name}.yaml"
    else:
        filename = f"{input_name}.{input_subtype}.yaml"
    return data_path / filename


class InvalidInputDefinition(Exception):
    pass


class MissingInputDefinition(InvalidInputDefinition):
    pass


def validate_input(input_params, required_fields):
    """
    Checks all entries in `input_params` against the definition in
    `required_fields`. The requirements may either be single value or function
    which validates the value provided and returns a deserialized version. To
    make conditional requirements the functions should accept the arguments
    `param_name` and `input_params` (the latter can the be queried to check the
    parameters the user has set).
    """

    def _check_field(f_name, f_option):
        """
        Validate an individual field with name `f_name` using the field
        definition given by `f_option` in the `input_params` dictionary.
        InvalidInputDefinition is raised if the value found in `input_params`
        is invalid. The cleaned value is returned, if None is returned then the
        value is optional and the user didn't provide a value (or they provided
        `None`)
        """
        # to enable implementation of conditional fields a special dictionary
        # with the fields `requires` and `choices` can be provided. `requires`
        # indicate other validations that must pass and choices are the valid
        # choices if they pass
        if type(f_option) == dict:
            if "requires" in f_option and "choices" in f_option:
                requirements = f_option["requires"]
                satisfied_requirements = {}
                for f_name_reqd, f_option_reqd in requirements.items():
                    # special flag for making a requirement that another
                    # parameter is set
                    if f_option_reqd == "__is_set__":
                        f_reqd_value = input_params.get(f_name_reqd)
                        f_value = input_params.get(f_name)
                        if f_reqd_value is not None:
                            satisfied_requirements[f_name_reqd] = f_value
                            continue
                        else:
                            raise InvalidInputDefinition(
                                f"For `{f_name}` == `{f_value}` the `{f_name_reqd}`"
                                " must be set"
                            )

                    try:
                        res = _check_field(f_name=f_name_reqd, f_option=f_option_reqd)
                        satisfied_requirements[f_name_reqd] = res
                    except InvalidInputDefinition:
                        break

                if len(satisfied_requirements) == len(requirements):
                    # all requirements satisfied, check the value provided
                    return _check_field(f_name=f_name, f_option=f_option["choices"])
                else:
                    # the requirements for this conditional parameter haven't
                    # been met and so the parameter shouldn't be set
                    if input_params.get(f_name) is not None:
                        raise InvalidInputDefinition(
                            f"`{f_name}` shouldn't be set unless the following "
                            f"requirements have been satisfied: `{requirements}`"
                        )
                    else:
                        return None
            else:
                raise NotImplementedError(
                    "A parameter validation provided as a dictionary must contain "
                    "the keys `requires` and `choices` to indicate the required "
                    "values of other fields and the valid choices when that "
                    "requirement is satisfied"
                )

        # Check for optional parameters (were a default value will
        # otherwise be used). This accomplished by providing `None` in the list
        # of options
        missing_allowed = type(f_option) in [list, tuple] and None in f_option
        if missing_allowed and f_name not in input_params:
            return None

        if f_name not in input_params:
            raise MissingInputDefinition("Missing `{}` field".format(f_name))

        if callable(f_option):
            # the callable will provide the correct conversion and also ensure
            # that the provided value has the correct type (otherwise we expect
            # the callable to raise an error)
            return f_option(input_params[f_name])

        # here we check whether the parameter was prescribed to have a specific
        # type and in that case whether the value provided has the correct type
        if type(f_option) == type and type(input_params[f_name]) != f_option:
            raise InvalidInputDefinition(
                "Field `{}` should have type {}, but has type `{}`".format(
                    f_name, f_option, type(input_params[f_name])
                )
            )

        # Check if the value provided is exactly the same as given by
        # `f_option`
        same_type = isinstance(input_params[f_name], type(f_option))
        same_value = input_params[f_name] == f_option
        if same_type and same_value:
            return f_option

        # choices a given as a list or tuple, so we call `_check_field`
        # recursively to see if one of the options fit
        if type(f_option) == list or type(f_option) == tuple:
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

        raise InvalidInputDefinition(f_name, f_option)

    checked_valid_fields = []
    for f_name, f_option in required_fields.items():
        if "|" in f_name:
            f_names = f_name.split("|")
            exceptions = []
            exceptions_missing_input = []
            valid_f_names = []
            for f_name in f_names:
                try:
                    new_val = _check_field(f_name, f_option)
                    if new_val is not None:
                        input_params[f_name] = new_val
                        valid_f_names.append(f_name)
                except MissingInputDefinition as e:
                    exceptions_missing_input.append(e)
                except Exception as e:
                    exceptions.append(e)

            # we allow for one missing part because `|` denotes an OR operation
            if len(exceptions) > 0 or len(exceptions_missing_input) > 1:
                raise InvalidInputDefinition(
                    "The following issues were found when trying to parse "
                    f"the values of {', '.join(f_names)} for the `{f_option}` parameters: "
                    f"{exceptions}"
                )
            else:
                checked_valid_fields += valid_f_names
        else:
            new_val = _check_field(f_name, f_option)
            if new_val is not None:
                input_params[f_name] = new_val
            elif f_name in input_params:
                del input_params[f_name]
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
