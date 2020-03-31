from pathlib import Path
import dateutil.parser
import datetime

import pytraj.produce.lagragian_trajectory


def test_create_stationary_trajectory():
    fn_out = "stationary_trajectory_test.nc"
    pytraj.produce.lagragian_trajectory.main(
        lat0=-10, lon0=40,
        t0=dateutil.parser.parse("2020-01-22T12:00Z"),
        t_max=dateutil.parser.parse("2020-01-24T12:00Z"),
        dt=datetime.timedelta(hours=4),
        trajectory_type='stationary',
        out_filename=fn_out,
    )

    p = Path(fn_out)
    assert p.exists()
    p.unlink()
