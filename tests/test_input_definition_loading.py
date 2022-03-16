import shutil
import tempfile
from pathlib import Path

from lagtraj.input_definitions import build_input_definition_path
from lagtraj.input_definitions.examples import EXAMPLES_ROOT_PATH as p_examples_root
from lagtraj.trajectory import load as traj_load


def test_loading_local_file_by_name():
    tempdir = tempfile.TemporaryDirectory()
    testdata_dir = Path(tempdir.name)

    input_name = "eurec4a_20200103_lag_testcase"
    input_type = "trajectory"

    # copy an example out from the lagtraj input examples
    p_example = build_input_definition_path(
        input_name=input_name,
        input_type=input_type,
        root_data_path=p_examples_root,
    )
    p_local = build_input_definition_path(
        input_name=input_name, input_type=input_type, root_data_path=testdata_dir
    )
    p_local.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(p_example, p_local)

    # and check we can load it - by name
    traj_load.load_definition(root_data_path=testdata_dir, name=input_name)
