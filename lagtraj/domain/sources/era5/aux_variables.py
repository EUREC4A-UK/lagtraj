import xarray as xr
import numpy as np

from .constants import rd, rv_over_rd_minus_one, cp, p_ref, rg, rlv, rls, omega
from ....utils.interpolation.methods import cos_transition
from ....utils.thermo import theta_l_detailed

rd_over_cp = rd / cp
p_ref_inv = 1.0 / p_ref
rlv_over_cp = rlv / cp
rls_over_cp = rls / cp


def calc_variable(ds, var, **kwargs):
    """
    Calculate an auxiliary variable `var` using ERA5 data in `ds`
    """
    if var == "theta":
        # Compute the potential temperature
        return xr.DataArray(
            ds.t * (ds.p_f * p_ref_inv) ** -rd_over_cp,
            attrs={"units": "K", "long_name": "potential temperature"},
        )
    elif var == "rho":
        # Computes density of moist air (without cloud water/ice and hydrometeors)
        return xr.DataArray(
            ds.p_f / (rd * ds.t * (1.0 + rv_over_rd_minus_one * ds.q)),
            attrs={"units": "kg m**-3", "long_name": "density"},
        )
    elif var == "w_pressure_corr":
        # Compute pressure velocity, corrected for diurnal cycle
        # with a correction that gradually vainishes between
        # pressures w_cutoff_start and w_cutoff_end
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
        # Compute the corrected vertical velocity in m/s
        return xr.DataArray(
            -ds.w_pressure_corr / (rg * ds.rho),
            attrs={"units": "m s**-1", "long_name": "Corrected vertical velocity"},
        )
    elif var == "t_l":
        # Compute the Liquid/ice water temperature
        return xr.DataArray(
            ds["t"] - rlv_over_cp * (ds["clwc"]) - rls_over_cp * (ds["ciwc"]),
            attrs={"units": "K", "long_name": "Liquid/ice water temperature"},
        )
    elif var == "q_t":
        # Compute the specific total water content without hydrometeors
        return xr.DataArray(
            ds["q"] + ds["clwc"] + ds["ciwc"],
            attrs={
                "long_name": "Specific total water content (no hydrometeors)",
                "units": "kg kg**-1",
            },
        )
    elif var == "q_t_hydromet":
        # Compute the specific total water content with hydrometeors
        return xr.DataArray(
            ds["q"] + ds["clwc"] + ds["ciwc"] + ds["crwc"] + ds["cswc"],
            attrs={
                "long_name": "Specific total water content (with hydrometeors)",
                "units": "kg kg**-1",
            },
        )
    elif var == "r_t":
        # Compute the total water mixing ratio without hydrometeors
        return xr.DataArray(
            ds["q_t"] / (1.0 - ds["q_t_hydromet"]),
            attrs={"long_name": "Total water mixing ratio", "units": "kg kg**-1"},
        )
    elif var == "r_v":
        # Compute the water vapour mixing ratio
        return xr.DataArray(
            ds["q"] / (1.0 - ds["q_t_hydromet"]),
            attrs={"long_name": "Water vapour mixing ratio", "units": "kg kg**-1"},
        )
    elif var == "r_l":
        # Compute the cloud liquid water mixing ratio
        return xr.DataArray(
            ds["clwc"] / (1.0 - ds["q_t_hydromet"]),
            attrs={
                "long_name": "Cloud liquid water mixing ratio",
                "units": "kg kg**-1",
            },
        )
    elif var == "r_i":
        # Compute the cloud ice vapour mixing ratio
        return xr.DataArray(
            ds["ciwc"] / (1.0 - ds["q_t_hydromet"]),
            attrs={"long_name": "Cloud ice mixing ratio", "units": "kg kg**-1"},
        )
    elif var == "theta_l":
        # Compute a relatively well conserved liquid water potential temperature
        return xr.DataArray(
            theta_l_detailed(ds["t"], ds["p_f"], ds["q_t"], ds["clwc"], ds["ciwc"],),
            attrs={"long_name": "Liquid water potential temperature", "units": "K"},
        )
    elif var == "u_g":
        f_cor = 2.0 * omega * np.sin(np.deg2rad(ds["lat"]))
        return xr.DataArray(
            -(1.0 / (f_cor * ds["rho_mean"])) * ds["dp_fdy"],
            attrs={"long_name": "U component of geostrophic wind", "units": "m s**-1"},
        )
    elif var == "v_g":
        f_cor = 2.0 * omega * np.sin(np.deg2rad(ds["lat"]))
        return xr.DataArray(
            (1.0 / (f_cor * ds["rho_mean"])) * ds["dp_fdx"],
            attrs={"long_name": "V component of geostrophic wind", "units": "m s**-1"},
        )
    else:
        raise NotImplementedError(f"Variable `{var}` not implemented")
