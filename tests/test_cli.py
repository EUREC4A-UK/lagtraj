from unittest.mock import patch
import importlib

import lagtraj.forcings.create
from lagtraj.forcings import build_forcing_data_path
import lagtraj.trajectory.create
from lagtraj.trajectory import build_data_path as build_trajectory_data_path

from test_forcing_profiles_extraction import AVAILABLE_CONVERSIONS


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
