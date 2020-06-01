from ..input_definitions import load
from . import INPUT_REQUIRED_FIELDS


def load_definition(data_path, forcing_name):
    return load.load_definition(
        data_path=data_path,
        input_name=forcing_name,
        input_type="forcing",
        required_fields=INPUT_REQUIRED_FIELDS,
    )
