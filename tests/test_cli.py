from unittest.mock import patch

import lagtraj.forcings.create
import lagtraj.trajectory.create

from test_forcing_profiles_extraction import AVAILABLE_CONVERSIONS


@patch("lagtraj.domain.download.download_complete")
def test_cli(download_complete_mocked, testdata_info):
    # during CI we won't have access to an ERA5 token so we need to pretend
    # that all data has been downloaded
    download_complete_mocked.return_value = True

    args = [
        testdata_info["trajectory_name"],
        "--data-path",
        str(testdata_info["testdata_path"]),
    ]
    lagtraj.trajectory.create.cli(args=args)

    for conversion in AVAILABLE_CONVERSIONS:
        args = [
            testdata_info["forcing_name"],
            "--data-path",
            str(testdata_info["testdata_path"]),
        ]

        if conversion is not None:
            args += ["--conversion", conversion]
        lagtraj.forcings.create.cli(args=args)
