import numpy as np


def test_check_domain_data(ds_domain_test):
    """
    Rudimentary tests to ensure that data used for testing is sensible
    """
    assert not np.isnan(ds_domain_test.isel(time=0).u.mean())
    assert not np.isnan(ds_domain_test.isel(time=-1).u.mean())
