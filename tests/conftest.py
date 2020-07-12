from pathlib import Path
import tarfile
import tempfile

import requests
import pytest

import lagtraj.domain.load

TESTDATA_URL = (
    "http://gws-access.ceda.ac.uk/public/eurec4auk/testdata/lagtraj.testdata.tar.gz"
)

# A testdata folder in this directory
tempdir = tempfile.TemporaryDirectory()
testdata_dir = Path(tempdir.name)


def download_testdata():
    fhtar = tempfile.NamedTemporaryFile(delete=False, suffix=".tar.gz")

    r = requests.get(TESTDATA_URL)
    fhtar.write(r.content)
    fhtar.close()

    testdata_dir.mkdir()
    tarfile.open(fhtar.name, "r:gz").extractall(testdata_dir)

    return


@pytest.fixture
def ds_domain_test(scope="session"):
    # Download testdata if it is not there yet
    if len(list(testdata_dir.glob("**/*.nc"))) == 0:
        print("Downloading testdata...")
        download_testdata()

    DOMAIN_NAME = "eurec4a_circle_eul"
    ds = lagtraj.domain.load.load_data(root_data_path=testdata_dir, name=DOMAIN_NAME)
    return ds
