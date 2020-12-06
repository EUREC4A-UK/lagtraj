import pytest

import lagtraj.input_definitions.examples
import lagtraj.input_definitions.load
from lagtraj import DATA_TYPE_PLURAL, DEFAULT_ROOT_DATA_PATH

from lagtraj.domain import INPUT_REQUIRED_FIELDS as DOMAIN_REQUIRED_FIELDS
from lagtraj.forcings import INPUT_REQUIRED_FIELDS as FORCING_REQUIRED_FIELDS
from lagtraj.trajectory import INPUT_REQUIRED_FIELDS as TRAJECTORY_REQUIRED_FIELDS
from lagtraj.conversion import INPUT_REQUIRED_FIELDS as CONVERSION_REQUIRED_FIELDS

INPUT_TYPES = {
    "forcing": FORCING_REQUIRED_FIELDS,
    "domain": DOMAIN_REQUIRED_FIELDS,
    "trajectory": TRAJECTORY_REQUIRED_FIELDS,
    "conversion": CONVERSION_REQUIRED_FIELDS,
}


def test_print_examples():
    for test_types in ["all", INPUT_TYPES.keys()]:
        lagtraj.input_definitions.examples.print_available(input_types=test_types)


def _get_examples():
    for input_type in INPUT_TYPES:

        def _visit_node(p):
            child_nodes = []
            for k, v in p.items():
                if len(v) == 0:
                    child_nodes.append(k)
                else:
                    child_nodes += [f"{k}/{cn}" for cn in _visit_node(v)]

            return child_nodes

        input_definition_examples = lagtraj.input_definitions.examples.get_available(
            input_types="all"
        )

        return _visit_node(input_definition_examples)


INPUT_DEFN_EXAMPLES = _get_examples()


@pytest.mark.parametrize("input_example", INPUT_DEFN_EXAMPLES)
def test_load_example(input_example):
    i = input_example.index("/")
    input_type_plural = input_example[:i]
    input_name = f"lagtraj://{input_example[i + 1 :]}"

    if input_type_plural == "forcings_conversion":
        # TODO: implement tests for the forcing conversions bundled with lagtraj
        return

    input_type = None
    for k, v in DATA_TYPE_PLURAL.items():
        if v == input_type_plural:
            input_type = k

    input_defn = lagtraj.input_definitions.load.load_definition(
        input_name=input_name,
        input_type=input_type,
        root_data_path=DEFAULT_ROOT_DATA_PATH,
        required_fields=INPUT_TYPES[input_type],
    )

    params = lagtraj.input_definitions.examples.attempt_read(
        input_name=input_name, input_type=input_type
    )

    # if any of the param values starts with `lagtraj://` we can assume this is
    # referring to another input-definition, so we should check here that it
    # exists
    for k, v in params.items():
        if isinstance(v, str) and k in INPUT_TYPES.keys():
            if not v.startswith("lagtraj://"):
                raise Exception("All input definitions which are included with "
                                "lagtraj should only refer to other inputs "
                                "included with lagtraj (i.e. should have the "
                                f"`lagtraj://` prefix). The `{input_name}` "
                                f"{input_type} input definition returns to the "
                                f"`{v}` {k}")
            else:
                try:
                    lagtraj.input_definitions.examples.attempt_read(
                        input_name=v, input_type=k
                    )
                except lagtraj.input_definitions.examples.LagtrajExampleDoesNotExist:
                    raise Exception(
                        f"The input-definition `lagtraj://{input_name}` "
                        f"refers to the `{v}` {k} input definition "
                        "which doesn't exist!"
                    )

    if input_defn["version"] == "unversioned":
        raise Exception(
            "All input definitions included with lagtraj should "
            "be versioned. Currently the input definition for "
            f"the `{input_name}` {input_type} is not"
        )
