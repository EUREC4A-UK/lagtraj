import xarray as xr

from .constants import rd, rv_over_rd_minus_one, cp, p_ref, rg
from ....utils.interpolation import cos_transition

rd_over_cp = rd / cp
p_ref_inv = 1.0 / p_ref


def calc_variable(ds, var, **kwargs):
    """
    Calculate an auxiliary variable `var` using ERA5 data in `ds`
    """
    if var == "theta":
        return xr.DataArray(
            ds.t * (ds.p_f * p_ref_inv) ** rd_over_cp,
            attrs={"units": "K", "long_name": "potential temperature"},
        )
    elif var == "rho":
        return xr.DataArray(
            ds.p_f / (rd * ds.t * (1.0 + rv_over_rd_minus_one * ds.q)),
            attrs={"units": "kg m**-3", "long_name": "density"},
        )
    elif var == "w_pressure_corr":
        if "w_cutoff_start" not in kwargs or "w_cutoff_end" not in kwargs:
            raise Exception(
                f"To calculate `{var}` the kwargs `w_cutoff_start`"
                " and `w_cutoff_end` are required"
            )
        # if not "p_f" in ds.data_vars:
        # p_f =
        vals = ds.w - ds.w[:, [-1], :, :].values * cos_transition(
            ds.p_f[:, :, :, :].values, kwargs["w_cutoff_start"], kwargs["w_cutoff_end"],
        )
        return xr.DataArray(
            vals,
            attrs={
                "units": ds.w.units,
                "long_name": "Corrected pressure vertical velocity",
            },
        )
    elif var == "w_corr":
        return xr.DataArray(
            -ds.w_pressure_corr / (rg * ds.rho),
            attrs={"units": "m s**-1", "long_name": "Corrected vertical velocity"},
        )
    else:
        raise NotImplementedError(f"Variable `{var}` not implemented")
