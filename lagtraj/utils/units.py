def fix_units(ds):
    """Changes units of data to make them compatible with the cf-checker"""
    units_dict = {
        "(0 - 1)": "1",
        "m of water equivalent": "m",
        "~": "1",
        "(0 - 1) s**-1": "s**-1",
        "(0 - 1) m**-1": "m**-1",
    }
    for variable in ds.variables:
        if hasattr(ds[variable], "units"):
            these_units = ds[variable].units
            if these_units in units_dict:
                ds[variable].attrs["units"] = units_dict[these_units]
