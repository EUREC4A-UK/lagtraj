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
    "conversion": CONVERSION_REQUIRED_FIELDS
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
    input_name = input_example[i + 1 :]

    input_type = None
    for k, v in DATA_TYPE_PLURAL.items():
        if v == input_type_plural:
            input_type = k

    lagtraj.input_definitions.load.load_definition(
        input_name=f"lagtraj://{input_name}",
        input_type=input_type,
        root_data_path=DEFAULT_ROOT_DATA_PATH,
        required_fields=INPUT_TYPES[input_type],
    )
