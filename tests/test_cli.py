import importlib
import re
from pathlib import Path
from unittest.mock import patch

import pytest
from conftest import TESTDATA_DIR, ensure_testdata_available
from test_forcing_profiles_extraction import AVAILABLE_CONVERSIONS

import lagtraj.forcings.create
import lagtraj.trajectory.create
from lagtraj.forcings import build_forcing_data_path
from lagtraj.trajectory import build_data_path as build_trajectory_data_path


@patch("lagtraj.domain.download.download_complete")
def test_cli(download_complete_mocked, testdata_info):
    # during CI we won't have access to an ERA5 token so we need to pretend
    # that all data has been downloaded. We have to reload the module so we
    # ensure the mocked function is called
    importlib.reload(lagtraj.trajectory.create)
    download_complete_mocked.return_value = True

    args = [
        testdata_info["trajectory_name"],
        "--data-path",
        str(testdata_info["testdata_path"]),
    ]
    lagtraj.trajectory.create.cli(args=args)

    # sort ensure `None` comes first
    conversions = sorted(AVAILABLE_CONVERSIONS, key=lambda x: (x is not None, x))

    # conversion creation where the non-converted forcing already exists,
    # beacuse `conversion=None` is run first
    for conversion in conversions:
        args = [
            testdata_info["forcing_name"],
            "--data-path",
            str(testdata_info["testdata_path"]),
        ]

        if conversion is not None:
            args += ["--conversion", conversion]
        lagtraj.forcings.create.cli(args=args)

    # cleanup, remove all saved forcings
    for conversion in conversions:
        p_forcing = build_forcing_data_path(
            root_data_path=testdata_info["testdata_path"],
            forcing_name=testdata_info["forcing_name"],
            conversion_name=conversion,
        )
        p_forcing.unlink()

    # conversion creation where the non-converted forcing doesn't exist
    for conversion in conversions:
        args = [
            testdata_info["forcing_name"],
            "--data-path",
            str(testdata_info["testdata_path"]),
        ]

        if conversion is not None:
            args += ["--conversion", conversion]
        lagtraj.forcings.create.cli(args=args)

        # delete saved forcing
        p_forcing = build_forcing_data_path(
            root_data_path=testdata_info["testdata_path"],
            forcing_name=testdata_info["forcing_name"],
            conversion_name=conversion,
        )
        p_forcing.unlink()

    # cleanup in case we run the CLI test again, the test will fail if we try
    # to overwrite an existing file
    p_traj = build_trajectory_data_path(
        root_data_path=testdata_info["testdata_path"],
        trajectory_name=testdata_info["trajectory_name"],
    )
    p_traj.unlink()


def _parse_readme_cli_commands():
    lines = open(Path(__file__).parent.parent / "README.md").read().splitlines()

    lines = [line.replace("$>", "").strip() for line in lines]

    cli_commands = list(filter(lambda line: "python -m lagtraj" in line, lines))

    return cli_commands


def _set_placeholder_args(cmd):
    """
    fill out placeholder args in README commands, replacing <required_arg> with
    sensible defaults and removing [--optional-arg], turning e.g.

        python -m lagtraj.forcing.create <forcing_name> [--conversion <conversion_name>]

    into

        python -m lagtraj.forcing.create lagtraj://eurec4a_20191209_12_lag [--conversion <conversion_name>]
    """
    # for examples where we've got placeholders we should replace these
    # placeholders with some reasonable defaults
    default_args = [
        ("<command>", "input_definitions.examples"),
        ("<domain_name>", "lagtraj://eurec4a_circle"),
        ("<trajectory_name>", "lagtraj://eurec4a_20200202_first_short"),
        ("<forcing_name>", "lagtraj://eurec4a_20200202_first_short"),
        ("<start_date>", "2020/02/02"),
        ("<end_date>", "2020/02/02"),
    ]
    for k, v in default_args:
        cmd = cmd.replace(k, v)

    for optional_arg in re.findall(r"\[.*\]", cmd):
        cmd = cmd.replace(optional_arg, "")

    return cmd


@pytest.fixture(scope="module", autouse=True)
def fetch_testdata():
    ensure_testdata_available()


@pytest.mark.parametrize("cli_command", _parse_readme_cli_commands())
def test_readme_cli_commands(cli_command):
    cli_command = _set_placeholder_args(cli_command)

    module_name, *args = cli_command.replace("python -m", "").strip().split(" ")

    try:
        module = importlib.import_module(module_name)
        cli_fn = module.cli
    except ImportError:
        raise NotImplementedError(
            f"Can't test for CLI command `{cli_command}` in README because it"
            " isn't clear what entrypoint function this command uses"
        )

    cmds_without_datapath = ["lagtraj.input_definitions.examples"]
    if module_name not in cmds_without_datapath:
        args += ["--data-path", str(TESTDATA_DIR.absolute())]

        # before we run the cli command we need to ensure that it only refers
        # to domain data that we know is included in the testdata dataset
        print(f"checking that we have testdata for command `{args}`")
        if not module.has_data_for_cli_command(args):
            raise Exception(
                f"The `{cli_command}` command in the README will use data "
                "which is not included in the testdata tar-ball"
            )

    print(f"running {module.__name__}.cli with args {args}")
    cli_fn(args)

    # quick hack to remove the output since otherwise forcing creation will
    # fail the second time it is run to create the same file
    [p.unlink() for p in (TESTDATA_DIR / "forcings").glob("*.nc")]
