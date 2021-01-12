import datetime


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


def round_time(dt=None, num_seconds=60):
    """Round a datetime object to any time lapse in seconds
    dt : datetime.datetime object, default now.
    num_seconds : Closest number of seconds to round to, default 1 minute.

    based off https://stackoverflow.com/a/10854034
    """
    if dt is None:
        dt = datetime.datetime.now()
    seconds = (dt.replace(tzinfo=None) - dt.min).seconds
    rounding = (seconds + num_seconds / 2) // num_seconds * num_seconds
    return dt + datetime.timedelta(0, rounding - seconds, -dt.microsecond)
