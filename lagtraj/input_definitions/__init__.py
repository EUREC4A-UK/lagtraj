from .. import build_data_path


def build_input_definition_path(root_data_path, input_name, input_type):
    data_path = build_data_path(root_data_path=root_data_path, data_type=input_type)

    return data_path / (input_name + ".yaml")


class InvalidInputDefinition(Exception):
    pass


def validate_input(input_params, required_fields):
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
                requirements = f_option['requires']
                satisfied_requirements = {}
                for f_name_reqd, f_option_reqd in requirements.items():
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
        elif f_name not in input_params:
            raise InvalidInputDefinition("Missing `{}` field".format(f_name))

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

    extra_fields = set(input_params.keys()).difference(checked_valid_fields)

    if len(extra_fields) > 0:
        raise InvalidInputDefinition(
            "Input definition has the following"
            " extra fields: {}".format(", ".join(extra_fields))
        )
