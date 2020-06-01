

INPUT_REQUIRED_FIELDS = dict(
    source=str,
    trajectory=str,
    domain=str,
    profiles_strategy=str,
    gradient_strategy=str,
    averaging_width=float,
    levels_strategy=str,
    levels_number=int,
    levels_dzmin=float,
    levels_top=float
)


class InvalidForcingDefinition(Exception):
    pass


def validate_input(domain_params):
    for f, f_type in INPUT_REQUIRED_FIELDS.items():
        if f not in domain_params:
            raise InvalidForcingDefinition("Missing `{}` field".format(f))
        elif type(domain_params[f]) != f_type:
            raise InvalidForcingDefinition(
                "Field `{}` should have type {}" "".format(f, f_type)
            )

    extra_fields = set(domain_params.keys()).difference(INPUT_REQUIRED_FIELDS.keys())

    if len(extra_fields) > 0:
        raise InvalidForcingDefinition(
            "Forcing definition has the following"
            " extra fields: {}".format(", ".join(extra_fields))
        )
