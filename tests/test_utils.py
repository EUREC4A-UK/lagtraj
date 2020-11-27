from lagtraj.utils.xarray import create_attributes_dictionary


def test_create_output_attributes():
    attrs_dict = dict(
        name="a kind of input",
        a_required_input=dict(
            name="name of required_input",
            param=42.0
        )
    )

    expected_attrs = dict(
        name=attrs_dict['name'],
        a_required_input_name=attrs_dict['a_required_input']['name'],
        a_required_input_param=attrs_dict['a_required_input']['param'],
    )

    actual_attrs = create_attributes_dictionary(attrs_dict)

    diffs = {}
    for k, v in expected_attrs.items():
        if not actual_attrs[k] == str(v):
            diffs[k] = (actual_attrs[k], v)

    assert len(diffs) == 0
