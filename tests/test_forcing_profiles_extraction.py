from pathlib import Path
import dateutil.parser
import datetime

import pytraj.produce.lagragian_trajectory
import pytraj.produce.forcing_profiles


def test_create_stationary_trajectory():
    fn_trajectory = "trajectory_test.nc"
    pytraj.produce.lagragian_trajectory.main(
        lat0=-10, lon0=40,
        t0=dateutil.parser.parse("2020-01-22T12:00"),
        t_max=dateutil.parser.parse("2020-01-24T12:00"),
        dt=datetime.timedelta(hours=4),
        trajectory_type='stationary',
        out_filename=fn_trajectory,
    )

    fn_forcing_profiles = "forcing_profiles.nc"
    pytraj.produce.forcing_profiles.main(
        fn_trajectory=fn_trajectory,
        source_data="dummy",
        out_filename=fn_forcing_profiles,
    )

    p_traj = Path(fn_trajectory)
    p_forcing = Path(fn_forcing_profiles)

    assert p_forcing.exists()
    p_traj.unlink()
    p_forcing.unlink()
