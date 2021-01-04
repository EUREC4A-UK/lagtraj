import pytest


from lagtraj.input_definitions import validate_input, InvalidInputDefinition


def test_single_value_parameter():
    required_fields = dict(
        a=float, b=int, c=str,
    )
    input_params = dict(
        a=2.0, b=42, c="foobar"
    )
    validate_input(input_params=input_params, required_fields=required_fields)


def test_choices_list():
    required_fields = dict(
        a=[None, float],
        b=["foo", "bar"]
    )
    validate_input(dict(a=2.0, b="foo"), required_fields=required_fields)
    validate_input(dict(b="bar"), required_fields=required_fields)

    with pytest.raises(InvalidInputDefinition):
        validate_input(dict(b="blah"), required_fields=required_fields)


def test_conditional():
    required_fields = dict(
        a=dict(requires=dict(b=["foo", "bar"]), choices=float),
        b=["foo", "bar", "blup"]
    )
    validate_input(dict(a=2.0, b="foo"), required_fields=required_fields)
    # a is only required when b is "foo" or "bar"
    validate_input(dict(b="blup"), required_fields=required_fields)

    with pytest.raises(InvalidInputDefinition):
        # when b == "foo" then a is required
        validate_input(dict(b="foo"), required_fields=required_fields)


def test_or():
    required_fields = {"a|b": float}
    validate_input(dict(a=2.0), required_fields=required_fields)
    validate_input(dict(b=2.0), required_fields=required_fields)
    validate_input(dict(a=42.0, b=2.0), required_fields=required_fields)

    with pytest.raises(InvalidInputDefinition):
        # one of a or b must be defined
        validate_input(dict(), required_fields=required_fields)
