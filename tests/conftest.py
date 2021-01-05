from pathlib import Path
import tarfile
import tempfile
import os

import requests
import pytest

import lagtraj.domain.load

TESTDATA_URL = "http://gws-access.jasmin.ac.uk/public/eurec4auk/testdata/lagtraj.testdata.v0.1.0.tar.gz"  # noqa

if os.environ.get("LAGTRAJ_TESTDATA_DIR", None):
    testdata_dir = Path(os.environ["LAGTRAJ_TESTDATA_DIR"])
else:
    tempdir = tempfile.TemporaryDirectory()
    testdata_dir = Path(tempdir.name)


def download_testdata():
    fhtar = tempfile.NamedTemporaryFile(delete=False, suffix=".tar.gz")

    r = requests.get(TESTDATA_URL)
    fhtar.write(r.content)
    fhtar.close()

    tarfile.open(fhtar.name, "r:gz").extractall(testdata_dir)


@pytest.fixture
def ds_domain_test(scope="session"):
    if not testdata_dir.exists():
        raise Exception(f"Couldn't find test-data directory {testdata_dir}")
    # Download testdata if it is not there yet
    if len(list(testdata_dir.glob("**/*.nc"))) == 0:
        print("Downloading testdata...")
        download_testdata()

    DOMAIN_NAME = "eurec4a_circle"
    ds = lagtraj.domain.load.load_data(root_data_path=testdata_dir, name=DOMAIN_NAME)
    return ds


@pytest.fixture
def ds_trajectory_linear(ds_domain_test):
    t0 = ds_domain_test.time.isel(time=-15)

    origin = lagtraj.trajectory.TrajectoryOrigin(
        lat=ds_domain_test.lat.mean(), lon=ds_domain_test.lon.mean(), datetime=t0,
    )

    da_times = ds_domain_test.time

    ds_traj = lagtraj.trajectory.create.create_trajectory(
        origin=origin, trajectory_type="linear", da_times=da_times, U=[0.0, -0.0]
    )

    return ds_traj
