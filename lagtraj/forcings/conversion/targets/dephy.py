"""This module contains:
    1. a list of the dephy variables as a dictionary
    2. the mapping between the internal lagtraj variables (era5) and dephy
    3. and a function to process the conversion
"""

import numpy as np
import xarray as xr


from ....domain.sources.era5.constants import cp, rlv
from ....utils.thermo import tref, qvsi, qvsl
from ....utils.interpolation.methods import (
    steffen_1d_no_ep_time,
    central_estimate,
    cos_transition,
)

EXPORT_KWARGS = dict(format="NETCDF3_CLASSIC")

dephy_variables = {
    "lat_traj": {"long_name": "Latitude", "units": "degrees_north"},
    "lon_traj": {"long_name": "Longitude", "units": "degrees_east"},
    "lat0_traj": {"long_name": "Initial latitude", "units": "degrees_north"},
    "lon0_traj": {"long_name": "Initial longitude", "units": "degrees_east"},
    "height": {"long_name": "Height above ground", "units": "m"},
    "pressure": {"long_name": "Pressure", "units": "Pa"},
    "u": {"long_name": "Zonal wind", "units": "m s-1"},
    "v": {"long_name": "Meridional wind", "units": "m s-1"},
    "temp": {"long_name": "Temperature", "units": "K"},
    "theta": {"long_name": "Potential temperature", "units": "K"},
    "thetal": {"long_name": "Liquid potential temperature", "units": "K"},
    "qv": {"long_name": "Specific humidity", "units": "kg kg-1"},
    "qt": {"long_name": "Total water content", "units": "kg kg-1"},
    "rv": {"long_name": "Water vapor mixing ratio", "units": "kg kg-1"},
    "rt": {"long_name": "Total water mixing ratio", "units": "kg kg-1"},
    "rl": {"long_name": "Liquid water mixing ratio", "units": "kg kg-1"},
    "ri": {"long_name": "Ice water mixing ratio", "units": "kg kg-1"},
    "ql": {"long_name": "Liquid water content", "units": "kg kg-1"},
    "qi": {"long_name": "Ice water content", "units": "kg kg-1"},
    "tke": {"long_name": "Turbulent kinetic energy", "units": "m2 s-2"},
    "time": {"long_name": "Forcing time"},
    "ps_forc": {"long_name": "Surface pressure for forcing", "units": "Pa"},
    "height_forc": {"long_name": "Height for forcing", "units": "m"},
    "pressure_forc": {"long_name": "Pressure for forcing", "units": "Pa"},
    "ug": {"long_name": "Geostrophic zonal wind", "units": "m s-1"},
    "vg": {"long_name": "Geostrophic meridional wind", "units": "m s-1"},
    "u_adv": {"long_name": "Zonal wind large-scale advection", "units": "m s-2"},
    "v_adv": {"long_name": "Meridional wind large-scale advection", "units": "m s-2"},
    "temp_adv": {"long_name": "Temperature large-scale advection", "units": "K s-1"},
    "theta_adv": {
        "long_name": "Potential temperature large-scale advection",
        "units": "K s-1",
    },
    "thetal_adv": {
        "long_name": "Liquid potential temperature large-scale advection",
        "units": "K s-1",
    },
    "temp_rad": {
        "long_name": "Radiative temperature large-scale tendency",
        "units": "K s-1",
    },
    "theta_rad": {
        "long_name": "Radiative potential temperature tendency",
        "units": "K s-1",
    },
    "thetal_rad": {
        "long_name": "Radiative liquid potential temperature tendency",
        "units": "K s-1",
    },
    "qv_adv": {
        "long_name": "Specific humidity large-scale advection",
        "units": "kg kg-1 s-1",
    },
    "qt_adv": {
        "long_name": "Total water content large-scale advection",
        "units": "kg kg-1 s-1",
    },
    "rv_adv": {
        "long_name": "Water vapor mixing ratio large-scale advection",
        "units": "kg kg-1 s-1",
    },
    "rt_adv": {
        "long_name": "Total water mixing ratio large-scale advection",
        "units": "kg kg-1 s-1",
    },
    "w": {"long_name": "Vertical velocity", "units": "m s-1"},
    "omega": {"long_name": "Pressure vertical velocity", "units": "Pa s-1"},
    "ts": {"long_name": "Surface temperature", "units": "K"},
    "ps": {"long_name": "Surface pressure", "units": "Pa"},
    "rh": {"long_name": "Relative humidity", "units": "%"},
    "temp_nudging": {"long_name": "Temperature for nudging", "units": "K"},
    "theta_nudging": {"long_name": "Potential temperature for nudging", "units": "K"},
    "thetal_nudging": {
        "long_name": "Liquid potential temperature for nudging",
        "units": "K",
    },
    "qv_nudging": {
        "long_name": "Specific humidity profile for nudging",
        "units": "kg kg-1",
    },
    "qt_nudging": {
        "long_name": "Total water content profile for nudging",
        "units": "kg kg-1",
    },
    "rv_nudging": {
        "long_name": "Water vapor mixing ratio profile for nudging",
        "units": "kg kg-1",
    },
    "rt_nudging": {
        "long_name": "Total water mixing ratio profile for nudging",
        "units": "kg kg-1",
    },
    "u_nudging": {"long_name": "Zonal wind profile for nudging", "units": "m s-1"},
    "v_nudging": {"long_name": "Meridional wind profile for nudging", "units": "m s-1"},
    "nudging_inv_u_traj": {
        "long_name": "Inverse nudging time profile for zonal wind",
        "units": "s-1",
    },
    "nudging_inv_v_traj": {
        "long_name": "Inverse nudging time profile for meridional wind",
        "units": "s-1",
    },
    "nudging_inv_temp_traj": {
        "long_name": "Inverse nudging time profile for temperature",
        "units": "s-1",
    },
    "nudging_inv_theta_traj": {
        "long_name": "Inverse nudging time profile for potential temperature",
        "units": "s-1",
    },
    "nudging_inv_thetal_traj": {
        "long_name": "Inverse nudging time profile for liquid potential temperature",
        "units": "s-1",
    },
    "nudging_inv_qv_traj": {
        "long_name": "Inverse nudging time profile for specific humidity",
        "units": "s-1",
    },
    "nudging_inv_qt_traj": {
        "long_name": "Inverse nudging time profile for total water content",
        "units": "s-1",
    },
    "nudging_inv_rv_traj": {
        "long_name": "Inverse nudging time profile for water vapor mixing ratio",
        "units": "s-1",
    },
    "nudging_inv_rt_traj": {
        "long_name": "Inverse nudging time profile for total water mixing ratio",
        "units": "s-1",
    },
    "sfc_sens_flx": {
        "long_name": "Surface sensible heat flux (positive upward)",
        "units": "W m-2",
    },
    "sfc_lat_flx": {
        "long_name": "Surface latent heat flux (positive upward)",
        "units": "W m-2",
    },
    "wpthetap": {
        "long_name": "Surface flux of potential temperature",
        "units": "K m s-1",
    },
    "wpqvp": {
        "long_name": "Surface flux of water vapor specific humidity",
        "units": "m s-1",
    },
    "wpqtp": {
        "long_name": "Surface flux of total water specific humidity",
        "units": "m s-1",
    },
    "wprvp": {
        "long_name": "Surface flux of water vapor mixing ratio",
        "units": "m s-1",
    },
    "wprtp": {
        "long_name": "Surface flux of total water mixing ratio",
        "units": "m s-1",
    },
    "z0_traj": {"units": "m", "long_name": "Roughness length for momentum"},
    "z0th_traj": {"units": "m", "long_name": "Roughness length for heat"},
    "ustar": {"units": "m s-1", "long_name": "Surface friction velocity"},
    "msnswrf": {
        "units": "W m-2",
        "long_name": "Mean surface net short-wave radiation flux",
    },
    "msnlwrf": {
        "units": "W m-2",
        "long_name": "Mean surface net long-wave radiation flux",
    },
    "mtnswrf": {
        "units": "W m-2",
        "long_name": "Mean top net short-wave radiation flux",
    },
    "mtnlwrf": {"units": "W m-2", "long_name": "Mean top net long-wave radiation flux"},
    "mtnswrfcs": {
        "units": "W m-2",
        "long_name": "Mean top net short-wave radiation flux, clear sky",
    },
    "mtnlwrfcs": {
        "units": "W m-2",
        "long_name": "Mean top net long-wave radiation flux, clear sky",
    },
    "msnswrfcs": {
        "units": "W m-2",
        "long_name": "Mean surface net short-wave radiation flux, clear sky",
    },
    "msnlwrfcs": {
        "units": "W m-2",
        "long_name": "Mean surface net long-wave radiation flux, clear sky",
    },
    "mtdwswrf": {
        "units": "W m-2",
        "long_name": "Mean top downward short-wave radiation flux",
    },
    "u_traj": {"units": "m s-1", "long_name": "Zonal wind speed of trajectory"},
    "v_traj": {"units": "m s-1", "long_name": "Meridional wind speed of trajectory"},
    "albedo_traj": {"units": "0-1", "long_name": "albedo"},
    "q_skin_traj": {"units": "m of water", "long_name": "skin reservoir content"},
}

# era5 units : dephy units
# we replace era5 units here
era5_to_dephy_units = {
    "m s**-1": "m s-1",
    "Pa s**-1": "Pa s-1",
    "m s**-1 s**-1": "m s-2",
    "metres": "m",
    "kg kg**-1": "kg kg-1",
    "K s**-1": "K s-1",
    "kg kg**-1 s**-1": "kg kg-1 s-1",
    "W m**-2": "W m-2",
    "1": "0-1",
}

# dephy variable : era5 variable
# we loop over dephy variables here
dephy_from_era5_initial_variables = {
    "lat0_traj": "lat",
    "lon0_traj": "lon",
    "height": "height_h_local",
    "pressure": "p_h_mean",
    "u": "u_mean",
    "v": "v_mean",
    "temp": "t_mean",
    "theta": "theta_mean",
    "thetal": "theta_l_mean",
    "qv": "q_mean",
    "qt": "q_t_mean",
    "rv": "r_v_mean",
    "rt": "r_t_mean",
    "rl": "r_l_mean",
    "ri": "r_i_mean",
    "ql": "clwc_mean",
    "qi": "ciwc_mean",
    "ps": "sp_mean",
}


# dephy variable : era5 variable
# we loop over dephy variables here
dephy_from_era5_forcing_variables = {
    "lat_traj": "lat",
    "lon_traj": "lon",
    "ps_forc": "sp_mean",
    "height_forc": "height_h_mean",
    "pressure_forc": "p_h_mean",
    "ug": "u_g",
    "vg": "v_g",
    "u_adv": "dudt_adv",
    "v_adv": "dvdt_adv",
    "temp_adv": "dtdt_adv",
    "theta_adv": "dthetadt_adv",
    "thetal_adv": "dtheta_ldt_adv",
    "qv_adv": "dqdt_adv",
    "qt_adv": "dq_tdt_adv",
    "rv_adv": "dr_vdt_adv",
    "rt_adv": "dr_tdt_adv",
    "w": "w_corr_mean",
    "omega": "w_pressure_corr_mean",
    "temp_nudging": "t_mean",
    "theta_nudging": "theta_mean",
    "thetal_nudging": "theta_l_mean",
    "qv_nudging": "q_mean",
    "qt_nudging": "q_t_mean",
    "rv_nudging": "r_v_mean",
    "rt_nudging": "r_t_mean",
    "u_nudging": "u_mean",
    "v_nudging": "v_mean",
    "ustar": "zust_mean",
    "z0_traj": "fsr_mean",
    "u_traj": "u_traj",
    "v_traj": "v_traj",
    "albedo_traj": "fal_mean",
}


def _rh_dephy(tt, pp, qt):
    """relative humidity, switching from liquid to ice at tref"""
    return xr.where(tt < tref, 100.0 * qt / qvsi(tt, pp), 100.0 * qt / qvsl(tt, pp))


def from_era5(ds_era5, da_levels, parameters, metadata):
    def init_field_dephy(field, variable):
        if np.ndim(field) == 1:
            return (
                ("t0", "lat", "lon"),
                field[:, None, None],
                dephy_variables[variable],
            )
        elif np.ndim(field) == 2:
            steffen_field = steffen_1d_no_ep_time(
                field, ds_era5["level"].values, dephy_level_array
            )
            return (
                ("t0", "lev", "lat", "lon"),
                steffen_field[:, :, None, None],
                dephy_variables[variable],
            )
        else:
            raise Exception("wrong dimension for dephy init field")

    def forcing_field_dephy(field, variable):
        if np.ndim(field) == 1:
            return (
                ("time", "lat", "lon"),
                field[:, None, None],
                dephy_variables[variable],
            )
        elif np.ndim(field) == 2:
            steffen_field = steffen_1d_no_ep_time(
                field, ds_era5["level"].values, dephy_level_array
            )
            return (
                ("time", "lev", "lat", "lon"),
                steffen_field[:, :, None, None],
                dephy_variables[variable],
            )
        else:
            raise Exception("wrong dimension for dephy forcing field")

    def unit_check(unit_guess, variable):
        if not unit_guess == dephy_variables[variable]["units"]:
            except_str = (
                "Incompatible units between ERA5 and dephy for variable "
                + variable
                + ". Please fix using the fix_era5_to_dephy_units dictionary:"
                + "ERA converted variable is "
                + unit_guess
                + ", dephy variable is "
                + dephy_variables[variable]["units"]
            )
            raise Exception(except_str)

    """Obtain a dephy input file from era5 variable set at high resolution"""
    dephy_level_array = da_levels.values
    dephy_level_coord = {
        "lev": ("lev", da_levels.values, {"long_name": "altitude", "units": "m"},),
    }
    dephy_t0_coord = {
        "t0": ("t0", [ds_era5.time.values[0]], {"long_name": "Initial time"})
    }
    dephy_time_coord = {
        "time": ("time", ds_era5.time.values, {"long_name": "Forcing time"})
    }
    dephy_lat_coord = {
        "lat": ("lat", [np.nan], {"long_name": "Latitude", "units": "degrees_north"})
    }
    dephy_lon_coord = {
        "lon": ("lon", [np.nan], {"long_name": "Longitude", "units": "degrees_east"})
    }
    ds_dephy = xr.Dataset(
        coords={
            **dephy_time_coord,
            **dephy_t0_coord,
            **dephy_level_coord,
            **dephy_lat_coord,
            **dephy_lon_coord,
        }
    )
    # Variables from dictionary
    # Including unit checks
    for variable in dephy_from_era5_initial_variables:
        era5_var = dephy_from_era5_initial_variables[variable]
        da_era5 = ds_era5[era5_var].isel(time=[0])
        # perform units check
        unit_guess = era5_to_dephy_units.get(da_era5.units, da_era5.units)
        unit_check(unit_guess, variable)
        ds_dephy[variable] = init_field_dephy(da_era5.values, variable)
    for variable in dephy_from_era5_forcing_variables:
        era5_var = dephy_from_era5_forcing_variables[variable]
        da_era5 = ds_era5[era5_var]
        # perform units check
        unit_guess = era5_to_dephy_units.get(da_era5.units, da_era5.units)
        unit_check(unit_guess, variable)
        ds_dephy[variable] = forcing_field_dephy(da_era5.values, variable)
    variables_to_centralise = {
        "msnswrf": "msnswrf_mean",
        "msnlwrf": "msnlwrf_mean",
        "mtnswrf": "mtnswrf_mean",
        "mtnlwrf": "mtnlwrf_mean",
        "mtnswrfcs": "mtnswrfcs_mean",
        "mtnlwrfcs": "mtnlwrfcs_mean",
        "msnswrfcs": "msnswrfcs_mean",
        "msnlwrfcs": "msnlwrfcs_mean",
        "mtdwswrf": "mtdwswrf_mean",
    }
    for variable in variables_to_centralise:
        this_central_estimate = central_estimate(
            ds_era5[variables_to_centralise[variable]].values
        )
        ds_dephy[variable] = forcing_field_dephy(this_central_estimate, variable)
    q_skin_field = ds_era5["src_mean"].values
    ds_dephy["q_skin_traj"] = forcing_field_dephy(q_skin_field, "q_skin_traj")
    # TKE (initial field) set to zero
    da_tke = 0.0 * ds_era5["u_mean"].isel(time=[0]) * ds_era5["u_mean"].isel(time=[0])
    ds_dephy["tke"] = init_field_dephy(da_tke.values, "tke")
    # Radiative tendencies, all-sky, combine SW and LW
    da_mtt = ds_era5["mttswr_mean"] + ds_era5["mttlwr_mean"]
    ds_dephy["temp_rad"] = forcing_field_dephy(da_mtt.values, "temp_rad")
    # In order to get theta/thetal tendencies, multiply by the exner function derived as theta/T
    da_mthetat = da_mtt * (ds_era5["theta_mean"] / ds_era5["t_mean"])
    ds_dephy["theta_rad"] = forcing_field_dephy(da_mthetat.values, "theta_rad")
    ds_dephy["thetal_rad"] = forcing_field_dephy(da_mthetat.values, "thetal_rad")
    # Heat roughness, derive from "flsr" variable
    z0th_traj = np.exp(ds_era5["flsr_mean"].values)
    ds_dephy["z0th_traj"] = forcing_field_dephy(z0th_traj, "z0th_traj")
    # Include same t_skin correction used for DALES, may need further work
    ts = ds_era5["stl1_mean"].values + 1.0
    ds_dephy["ts"] = forcing_field_dephy(ts, "ts")
    # Surface fluxes: obtain from time mean in ERA data, change sign for dephy!
    sfc_sens_flx = -central_estimate(ds_era5["msshf_mean"].values)
    ds_dephy["sfc_sens_flx"] = forcing_field_dephy(sfc_sens_flx, "sfc_sens_flx")
    sfc_lat_flx = -central_estimate(ds_era5["mslhf_mean"].values)
    ds_dephy["sfc_lat_flx"] = forcing_field_dephy(sfc_lat_flx, "sfc_lat_flx")
    wpthetap = (
        (sfc_sens_flx / (cp * ds_era5["rho_mean"].sel(level=0.0)))
        * (ds_era5["theta_mean"].sel(level=0.0) / ds_era5["t_mean"].sel(level=0.0))
    ).values
    ds_dephy["wpthetap"] = forcing_field_dephy(wpthetap, "wpthetap")
    wpqvp = (sfc_lat_flx / (rlv * ds_era5["rho_mean"].sel(level=0.0))).values
    ds_dephy["wpqvp"] = forcing_field_dephy(wpqvp, "wpqvp")
    wpqtp = wpqvp
    ds_dephy["wpqtp"] = forcing_field_dephy(wpqtp, "wpqtp")
    # Ratio of fluxes (mixing ratio vs. specific humidity) is same for all fluxes
    moisture_ratio = (
        ds_era5["r_t_local"].sel(level=0.0) / ds_era5["q_t_local"].sel(level=0.0)
    ).values
    wprvp = wpqvp * moisture_ratio
    ds_dephy["wprvp"] = forcing_field_dephy(wprvp, "wprvp")
    wprtp = wpqtp * moisture_ratio
    ds_dephy["wprtp"] = forcing_field_dephy(wprtp, "wprtp")
    rh = _rh_dephy(ds_dephy["temp"], ds_dephy["pressure"], ds_dephy["qt"])
    ds_dephy["rh"] = init_field_dephy(rh.values[:, :, 0, 0], "rh")

    def nudging_inv_time_prof(nudging_parameters, variable):
        height_array = ds_dephy["height_forc"].values
        nudging_method = nudging_parameters.method
        nudging_time = nudging_parameters.time
        nudging_height = nudging_parameters.height
        nudging_transition = nudging_parameters.transition
        if nudging_method is None:
            inv_time_array = 0.0 * height_array
        elif nudging_method == "cos":
            div_factor = cos_transition(
                height_array,
                nudging_height + 0.5 * nudging_transition,
                nudging_height - 0.5 * nudging_transition,
            )
            inv_time_array = div_factor / nudging_time
        else:
            raise Exception(
                "Method for calculating inverse nudging time profile undefined"
            )
        return (
            ("time", "lev", "lat", "lon"),
            inv_time_array,
            dephy_variables[variable],
        )

    nudging_parameters_momentum_traj = parameters.nudging_parameters_momentum_traj
    nudging_parameters_scalar_traj = parameters.nudging_parameters_scalar_traj
    ds_dephy["nudging_inv_u_traj"] = nudging_inv_time_prof(
        nudging_parameters_momentum_traj, "nudging_inv_u_traj"
    )
    ds_dephy["nudging_inv_v_traj"] = nudging_inv_time_prof(
        nudging_parameters_momentum_traj, "nudging_inv_v_traj"
    )
    ds_dephy["nudging_inv_temp_traj"] = nudging_inv_time_prof(
        nudging_parameters_scalar_traj, "nudging_inv_temp_traj"
    )
    ds_dephy["nudging_inv_theta_traj"] = nudging_inv_time_prof(
        nudging_parameters_scalar_traj, "nudging_inv_theta_traj"
    )
    ds_dephy["nudging_inv_thetal_traj"] = nudging_inv_time_prof(
        nudging_parameters_scalar_traj, "nudging_inv_thetal_traj"
    )
    ds_dephy["nudging_inv_qv_traj"] = nudging_inv_time_prof(
        nudging_parameters_scalar_traj, "nudging_inv_qv_traj"
    )
    ds_dephy["nudging_inv_qt_traj"] = nudging_inv_time_prof(
        nudging_parameters_scalar_traj, "nudging_inv_qt_traj"
    )
    ds_dephy["nudging_inv_rv_traj"] = nudging_inv_time_prof(
        nudging_parameters_scalar_traj, "nudging_inv_rv_traj"
    )
    ds_dephy["nudging_inv_rt_traj"] = nudging_inv_time_prof(
        nudging_parameters_scalar_traj, "nudging_inv_rt_traj"
    )

    # Final checks: are all variables present?
    for var in dephy_variables:
        if var not in ds_dephy:
            print(var + " is missing in the dephy formatted output")
    # Needs improvement

    dephy_dictionary = {
        "Conventions": "CF-1.0",
        "comment": metadata.comment,
        "reference": metadata.reference,
        "author": metadata.author,
        "modifications": metadata.modifications,
        "case": metadata.campaign + " " + metadata.case,
        "script": "https://github.com/EUREC4A-UK/lagtraj",
        "startDate": ds_dephy["time"][0].values.astype("str"),
        "endDate": ds_dephy["time"][-1].values.astype("str"),
        "adv_temp": parameters.adv_temp,
        "adv_theta": parameters.adv_theta,
        "adv_thetal": parameters.adv_thetal,
        "adv_qv": parameters.adv_qv,
        "adv_qt": parameters.adv_qt,
        "adv_rv": parameters.adv_rv,
        "adv_rt": parameters.adv_rt,
        "rad_temp": parameters.rad_temp,
        "rad_theta": parameters.rad_theta,
        "rad_thetal": parameters.rad_thetal,
        "forc_omega": parameters.forc_omega,
        "forc_w": parameters.forc_w,
        "forc_geo": parameters.forc_geo,
        "nudging_u": parameters.nudging_u,
        "nudging_v": parameters.nudging_v,
        "nudging_temp": parameters.nudging_temp,
        "nudging_theta": parameters.nudging_theta,
        "nudging_thetal": parameters.nudging_thetal,
        "nudging_qv": parameters.nudging_qv,
        "nudging_qt": parameters.nudging_qt,
        "nudging_rv": parameters.nudging_rv,
        "nudging_rt": parameters.nudging_rt,
        "z_nudging_u": np.nan,
        "z_nudging_v": np.nan,
        "z_nudging_temp": np.nan,
        "z_nudging_theta": np.nan,
        "z_nudging_thetal": np.nan,
        "z_nudging_qv": np.nan,
        "z_nudging_qt": np.nan,
        "z_nudging_rv": np.nan,
        "z_nudging_rt": np.nan,
        "p_nudging_u": np.nan,
        "p_nudging_v": np.nan,
        "p_nudging_temp": np.nan,
        "p_nudging_theta": np.nan,
        "p_nudging_thetal": np.nan,
        "p_nudging_qv": np.nan,
        "p_nudging_qt": np.nan,
        "p_nudging_rv": np.nan,
        "p_nudging_rt": np.nan,
        "zorog": 0.0,
        "z0": np.nan,
        "surfaceType": parameters.surfaceType,
        "surfaceForcing": parameters.surfaceForcing,
        "surfaceForcingWind": parameters.surfaceForcingWind,
    }
    ds_dephy.attrs.update(dephy_dictionary)
    return ds_dephy
