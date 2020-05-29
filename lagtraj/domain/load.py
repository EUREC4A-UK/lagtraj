from pathlib import Path

import xarray as xr


def load(data_path, name):
    domain_data_path = Path(data_path)/name

    ds = xr.open_mfdataset(domain_data_path/"*.nc")

    return ds
