"""
For nudging parameters we have quite a few sets of parameters and they can be
used to set the parameters for nudging momentum and scalar forcing profiles
separately (for example `nudging_timescale_scalars` and
`nudging_timescale_momentum`). To avoid typos while enumerating all those
options (and the cases for each nudging method were each parameter is required)
this code takes care of generating the valid combinations.
"""


AVAILABLE_NUDGING_METHODS = [
    "fixed_height",
    "runtime_inversion_height",
    "constant",
]

NUDGING_PARAMS_VALID_VALUES = dict(
    height=float, layer_thickness=float, timescale=int, shape=["cos"]
)


def _construct_nudging_required_fields():
    """
    parameters for nudging scalars and momentum are set separately, so we add
    them here in turn to the set of required parameters
    """
    fields = {}
    for nudged_variable_groups in ["scalars", "momentum"]:
        fields[
            f"nudging_method_{nudged_variable_groups}"
        ] = AVAILABLE_NUDGING_METHODS + [
            None,
        ]

        for param, param_choices in NUDGING_PARAMS_VALID_VALUES.items():
            param_fullname = f"nudging_{param}_{nudged_variable_groups}"

            if param == "height":
                required_for_nudging_method = ["fixed_height"]
            elif param == "timescale":
                required_for_nudging_method = AVAILABLE_NUDGING_METHODS
            elif param in ["layer_thickness", "shape"]:
                required_for_nudging_method = [
                    "fixed_height",
                    "runtime_inversion_height",
                ]
            else:
                raise NotImplementedError(param)

            fields[param_fullname] = dict(
                requires={
                    f"nudging_method_{nudged_variable_groups}": required_for_nudging_method,
                },
                choices=param_choices,
            )

    return fields


NUDGING_REQUIRED_FIELDS = _construct_nudging_required_fields()
