import numpy as np


def ds_time_to_seconds(ds):
    """Use seconds rather than hours as time unit in data"""
    attrs = {
        "long_name": "time",
        "units": "seconds since " + ds.time[0].values.astype("str"),
        "calendar": "proleptic_gregorian",
    }
    # Not sure why a float is needed
    ds["time"] = (
        "time",
        ((ds["time"] - ds.time[0]) / np.timedelta64(1, "s")).values,
        attrs,
    )
