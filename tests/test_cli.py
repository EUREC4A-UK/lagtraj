import lagtraj.forcings.create
import lagtraj.trajectory.create

from test_forcing_profiles_extraction import AVAILABLE_CONVERSIONS


def test_cli(testdata_info):
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
