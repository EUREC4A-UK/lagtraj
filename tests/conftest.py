from pathlib import Path
import tarfile
import tempfile
import os

import requests
import pytest

import lagtraj.domain.load

TESTDATA_URL = (
    "http://gws-access.ceda.ac.uk/public/eurec4auk/testdata/lagtraj.testdata.tar.gz"
)

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

    DOMAIN_NAME = "eurec4a_circle_eul"
    ds = lagtraj.domain.load.load_data(root_data_path=testdata_dir, name=DOMAIN_NAME)
    return ds
