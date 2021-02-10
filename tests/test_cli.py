import lagtraj.forcings.create
import lagtraj.trajectory.create

from test_forcing_profiles_extraction import AVAILABLE_CONVERSIONS


def test_cli(testdata_info):
    lagtraj.trajectory.create.cli(
        data_path=testdata_info["testdata_path"],
        trajectory_name=testdata_info["trajectory_name"],
    )

    for conversion in AVAILABLE_CONVERSIONS:
        lagtraj.forcings.create.cli(
            data_path=testdata_info["testdata_path"],
            forcing_name=testdata_info["forcing_name"],
            conversion_name=conversion,
        )
