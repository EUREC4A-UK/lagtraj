"""
Ensure that output of lagtraj doesn't produce errors in ncview by setting the
correct encoding for output and running ncview to check its output
"""
import tempfile
import signal
import psutil
import subprocess
import os


NCVIEW_COMPATIBLE_NETCDF_ENCODING = dict(
    time=dict(dtype="float64", units="seconds since 1970-01-01 00:00:00")
)


def build_valid_encoding(ds):
    encoding = {}
    for v in list(ds.data_vars) + list(ds.coords):
        try:
            ds[
                v
            ].dt  # try accessing the datetime accessor, only datetime types have this
            encoding[v] = dict(NCVIEW_COMPATIBLE_NETCDF_ENCODING["time"])
        except TypeError:
            pass
    return encoding


PROBLEM_STRINGS = [
    "ncview: netcdf_dim_value: unknown data type",
]


def check_for_ncview_warnings(ds):
    """
    Check for warnings when opening `ds` written to at netCDF file
    """
    fhnc = tempfile.NamedTemporaryFile(delete=False, suffix=".nc")
    encoding = build_valid_encoding(ds=ds)
    ds.to_netcdf(fhnc, encoding=encoding)
    fhnc.close()

    nc_filename = fhnc.name

    call = f"ncview {nc_filename}"
    p = subprocess.Popen(
        call.split(" "),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    def get_cpu_usage():
        # watch the ncview process for a second and see how much cpu it's using
        return psutil.Process(p.pid).cpu_percent(interval=1.0)

    # wait until ncview has finished loading
    while get_cpu_usage() > 0.1:
        continue

    # then kill it
    p.send_signal(signal.SIGINT)

    # appears ncview writes all its output to stderr...
    _, stderr = p.communicate()

    if any([s in stderr for s in PROBLEM_STRINGS]):
        raise Exception(f"ncview is not happy with the provided dataset: {stderr}")
