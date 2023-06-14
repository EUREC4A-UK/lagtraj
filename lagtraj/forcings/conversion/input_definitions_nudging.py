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
    False,
]

NUDGING_PARAMS_VALID_VALUES = dict(
    above_height=float,
    transition_thickness=float,
    timescale=int,
    transition_shape=["cos"],
)


def _construct_nudging_required_parameters():
    """
    parameters for nudging scalars and momentum are set separately, so we add
    them here in turn to the set of required parameters
    """
    nudging_parameters = {}
    for nudged_variable_groups in ["scalars", "momentum"]:
        nudging_parameters[
            f"nudging_method_{nudged_variable_groups}"
        ] = AVAILABLE_NUDGING_METHODS + [
            None,
        ]

        for param, param_choices in NUDGING_PARAMS_VALID_VALUES.items():
            param_fullname = f"nudging_{param}_{nudged_variable_groups}"

            if param == "above_height":
                required_for_nudging_method = [
                    "fixed_height",
                    "runtime_inversion_height",
                ]
            elif param == "timescale":
                required_for_nudging_method = [
                    "fixed_height",
                    "runtime_inversion_height",
                    "constant",
                ]
            elif param in ["transition_thickness", "transition_shape"]:
                required_for_nudging_method = [
                    "fixed_height",
                    "runtime_inversion_height",
                ]
            else:
                raise NotImplementedError(param)

            nudging_parameters[param_fullname] = dict(
                requires={
                    f"nudging_method_{nudged_variable_groups}": required_for_nudging_method,
                },
                choices=param_choices,
            )
    return nudging_parameters


NUDGING_REQUIRED_FIELDS = _construct_nudging_required_parameters()
