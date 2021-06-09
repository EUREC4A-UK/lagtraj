from unittest.mock import patch
import importlib
from pathlib import Path
import pytest
import subprocess
import re

import lagtraj.forcings.create
from lagtraj.forcings import build_forcing_data_path
import lagtraj.trajectory.create
from lagtraj.trajectory import build_data_path as build_trajectory_data_path

from test_forcing_profiles_extraction import AVAILABLE_CONVERSIONS


def _execute(cmd):
    # https://stackoverflow.com/a/4417735
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()

    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


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

    cli_commands = list(filter(lambda l: "python -m" in l, lines))

    # skip any commands with contain `[]` or `<>` pairs
    def has_optional_param(cmd):
        return len(re.findall(r"<.*>", cmd) + re.findall(r"\[.*\]", cmd)) > 0

    cli_commands = list(filter(lambda cmd: not has_optional_param(cmd), cli_commands))

    return cli_commands


@pytest.mark.parametrize("cli_command", _parse_readme_cli_commands())
@patch("lagtraj.domain.download.download_complete")
def test_readme_cli_commands(download_complete_mocked, cli_command):
    # during CI we won't have access to an ERA5 token so we need to pretend
    # that all data has been downloaded. We have to reload the module so we
    # ensure the mocked function is called
    importlib.reload(lagtraj.trajectory.create)
    download_complete_mocked.return_value = True

    for output in _execute(cli_command.split()):
        print((output.strip()))
