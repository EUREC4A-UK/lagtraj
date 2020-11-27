""" Thermodynamics routine using definitions based on the work by
Romps and Kuang (2010), but with specific humidities rather than
mixing ratios.
"""
import xarray as xr
import numpy as np

cpd = 1005.7  # specific heat at constant pressure (dry air).
cpv = 1870.0  # specific heat at constant pressure (vapour).
cpl = 4190.0  # specific heat at constant pressure (liquid).
cpi = 2106.0  # specific heat at constant pressure (ice).
srd = 6775.0  # specific entropy at triple point (dry air).
srv = 10320.0  # specific entropy at triple point (vapour).
srl = 3517.0  # specific entropy at triple point (liquid).
sri = 2296.0  # specific entropy at triple point (ice).
hrd = 530000.0  # specific enthalpy at triple point (dry air).
hrv = 3132927.0  # specific enthalpy at triple point (vapour).
hrl = 632000.0  # specific enthalpy at triple point (liquid).
hri = 298484.0  # specific enthalpy at triple point (ice).
rd = 287.04  # gas constant for dry air.
rv = 461.5  # gas constant for water vapor.
tref = 273.15  # reference temperature.
pref = 100000.0  # reference pressure.
p_triple_point = 610.7  # triple point pressure


def theta_l_detailed(tt, pp, qt, ql, qi):
    """theta_l: becomes theta for a dry parcel
    The default calculation is used whenever there is significant qv
    In the absence of qv but presence of ql (usually very high in the atmosphere,
    a correction is made to prevent division by zero."""
    theta_l = xr.where(
        ql + qi < 0.999 * qt,
        (
            tt
            * (tt / tref)
            ** (
                qt * (cpv - cpd) / cpd + ql * (cpl - cpv) / cpd + qi * (cpi - cpv) / cpd
            )
            * (pp / (pref * (1 + ((qt - ql - qi) * rv) / ((1 - qt) * rd))))
            ** (-(1 - qt) * rd / cpd)
            * (pp / (pref * (1 + ((1 - qt) * rd) / ((qt - ql - qi) * rv))))
            ** (-(qt - ql - qi) * rv / cpd)
            * np.exp(-ql * (srv - srl) / cpd - qi * (srv - sri) / cpd)
            * (1.0 / (1.0 + (qt * rv) / ((1.0 - qt) * rd))) ** ((1.0 - qt) * rd / cpd)
            * (1.0 / (1.0 + ((1.0 - qt) * rd) / (qt * rv))) ** ((qt) * rv / cpd)
        ),
        (
            tt
            * (tt / tref)
            ** (
                qt * (cpv - cpd) / cpd + ql * (cpl - cpv) / cpd + qi * (cpi - cpv) / cpd
            )
            * (pp / (pref * (1 + ((qt - ql - qi) * rv) / ((1 - qt) * rd))))
            ** (-(1 - qt) * rd / cpd)
            * np.exp(-ql * (srv - srl) / cpd - qi * (srv - sri) / cpd)
            * (1.0 / (1.0 + (qt * rv) / ((1.0 - qt) * rd))) ** ((1.0 - qt) * rd / cpd)
            * (1.0 / (1.0 + ((1.0 - qt) * rd) / (qt * rv))) ** ((qt) * rv / cpd)
        ),
    )
    return theta_l


def esatl(tt):
    """saturation pressure over the liquid phase"""
    return p_triple_point * np.exp(
        (1.0 / tref - 1.0 / tt) * ((hrv - hrl) - (cpv - cpl) * tref) / rv
        + ((cpv - cpl) / rv) * np.log(tt / tref)
    )


def esati(tt):
    """saturation pressure over the ice phase"""
    return p_triple_point * np.exp(
        (1.0 / tref - 1.0 / tt) * ((hrv - hri) - (cpv - cpi) * tref) / rv
        + ((cpv - cpi) / rv) * np.log(tt / tref)
    )


def qvsl(tt, pp):
    """saturation specific humidity over the liquid phase"""
    return rd / rv * esatl(tt) / (pp - (1.0 - rd / rv) * esatl(tt))


def qvsi(tt, pp):
    """saturation specific humidity over the ice phase"""
    return rd / rv * esati(tt) / (pp - (1.0 - rd / rv) * esati(tt))
