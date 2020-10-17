"""Module that enables conversions from era5 format to other formats
"""

from pathlib import Path
import datetime
import xarray as xr
import numpy as np

from .. import DEFAULT_ROOT_DATA_PATH
from ..forcings.load import load_data as load_forcing_data
from . import load, build_conversion_data_path
from ..utils import optional_debugging, validation
from ..utils.thermo import rh_dephy
from ..utils.interpolation import (
    steffen_1d_no_ep_time,
    central_estimate,
    cos_transition,
)
from .utils.levels import make_levels
from .kpt import kpt_variables
from .dephy import dephy_variables
from ..domain.sources.era5.constants import rg, cp, rlv


# TODO:
# - Check forcing has the right format?
# - Fix missing variables?
# Use ..utils.xarray import create_attributes_dictionary
# Add nudging


def kpt_from_era5(ds_era5, da_levels, parameters, metadata):
    """Obtain a kpt input file from era5 variable set at high resolution"""
    # Put full levels midway between half-levels, I think this is consistent with DALES
    # Reverse order of data, to confirm to other kpt input
    kpt_half_level_array = da_levels.values
    kpt_full_level_array = 0.5 * (
        kpt_half_level_array[:-1] + kpt_half_level_array[1:]
    )
    kpt_half_level_coord = {
        "nlevp1": (
            "nlevp1",
            (np.arange(len(kpt_half_level_array)) + 1.0)[::-1],
            {"long_name": "model half levels"},
        )
    }
    kpt_full_level_coord = {
        "nlev": (
            "nlev",
            (np.arange(len(kpt_full_level_array)) + 1.0)[::-1],
            {"long_name": "model full levels"},
        )
    }
    nDS_coord = {"nDS": ("nDS", [0], {},)}
    kpt_soil_coord = {
        "nlevs": ("nlevs", np.arange(4) + 1.0, {"long_name": "soil levels"})
    }
    ds_kpt = xr.Dataset(
        coords={
            "time": ds_era5.time,
            **kpt_full_level_coord,
            **kpt_half_level_coord,
            **kpt_soil_coord,
            **nDS_coord,
        }
    )
    # Variables from dictionary
    # Including unit checks
    for variable in kpt_from_era5_variables:
        era5_var = kpt_from_era5_variables[variable]
        da_era5 = ds_era5[era5_var]
        # perform units check
        unit_guess = era5_to_kpt_units.get(da_era5.units, da_era5.units)
        if not unit_guess == kpt_variables[variable]["units"]:
            except_str = (
                "Incompatible units between ERA5 and kpt for variable "
                + variable
                + ". Please fix using the fix_era5_to_kpt_units "
                + "dictionary: ERA converted variable is "
                + unit_guess
                + ", kpt variable is "
                + kpt_variables[variable]["units"]
            )
            raise Exception(except_str)
        # single level variable
        if np.ndim(da_era5.values) == 1:
            ds_kpt[variable] = (("time"), da_era5, kpt_variables[variable])
        # half level variable
        elif variable in ["zh", "presh"]:
            da_era5_on_half_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, kpt_half_level_array
            )
            ds_kpt[variable] = (
                ("time", "nlevp1"),
                da_era5_on_half_levels,
                kpt_variables[variable],
            )
        # full level variable
        else:
            da_era5_on_full_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, kpt_full_level_array
            )
            ds_kpt[variable] = (
                ("time", "nlev"),
                da_era5_on_full_levels,
                kpt_variables[variable],
            )
    # Simple unit fix fails for these variables
    # So these are added manually after checking
    # that units are compatible
    variables_to_manually_add = {
        "high_veg_type": "tvh_mean",
        "low_veg_type": "tvl_mean",
        "high_veg_lai": "lai_hv_mean",
        "low_veg_lai": "lai_lv_mean",
        "slor": "slor_mean",
        "q_skin": "src_mean",
        "snow": "sd_mean",
        "lsm": "lsm_mean",
    }
    for variable in variables_to_manually_add:
        ds_kpt[variable] = ds_era5[variables_to_manually_add[variable]]
        ds_kpt[variable].attrs.update(**kpt_variables[variable])
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
        ds_kpt[variable] = (
            ("time"),
            this_central_estimate,
            kpt_variables[variable],
        )
    # Soil moisture: combine levels
    swvl1 = ds_era5["swvl1_mean"].values
    swvl2 = ds_era5["swvl2_mean"].values
    swvl3 = ds_era5["swvl3_mean"].values
    swvl4 = ds_era5["swvl4_mean"].values
    q_soil = np.stack((swvl1, swvl2, swvl3, swvl4), axis=-1)
    ds_kpt["q_soil"] = (
        ("time", "nlevs"),
        q_soil,
        kpt_variables["q_soil"],
    )
    # Soil temperature: combine levels
    stl1 = ds_era5["stl1_mean"].values
    stl2 = ds_era5["stl2_mean"].values
    stl3 = ds_era5["stl3_mean"].values
    stl4 = ds_era5["stl4_mean"].values
    t_soil = np.stack((stl1, stl2, stl3, stl4), axis=-1)
    ds_kpt["t_soil"] = (
        ("time", "nlevs"),
        t_soil,
        kpt_variables["t_soil"],
    )
    # Soil thickness: combine levels
    h_soil = np.array([0.07, 0.21, 0.72, 1.89])
    ds_kpt["h_soil"] = (
        ("nlevs"),
        h_soil,
        kpt_variables["h_soil"],
    )
    # Orography: derive from surface geopotential
    ds_kpt["orog"] = ds_era5["z_mean"] / rg
    ds_kpt["orog"].attrs.update(**kpt_variables["orog"])
    # Heat roughness, derive from "flsr" variable
    ds_kpt["heat_rough"] = np.exp(ds_era5["flsr_mean"])
    ds_kpt["heat_rough"].attrs.update(**kpt_variables["heat_rough"])
    # Apply correction to t_skin (see output files)
    ds_kpt["t_skin"] = ds_era5["stl1_mean"] + 1.0
    ds_kpt["t_skin"].attrs.update(**kpt_variables["t_skin"])
    # Surface fluxes: obtain from time mean in ERA data, do not change sign
    sfc_sens_flx = central_estimate(ds_era5["msshf_mean"].values)
    ds_kpt["sfc_sens_flx"] = (
        ("time"),
        sfc_sens_flx,
        kpt_variables["sfc_sens_flx"],
    )
    sfc_lat_flx = central_estimate(ds_era5["mslhf_mean"].values)
    ds_kpt["sfc_lat_flx"] = (
        ("time"),
        sfc_lat_flx,
        kpt_variables["sfc_lat_flx"],
    )
    # Final checks: are all variables present?
    ds_kpt["time_traj"] = (
        ds_era5["time"] - np.datetime64("1970-01-01T00:00")
    ) / np.timedelta64(1, "s")
    ds_kpt["time_traj"].attrs.update(**kpt_variables["time_traj"])
    ds_kpt["DS"] = (("nDS"), ["Trajectory origin"], kpt_variables["DS"])
    ds_kpt["timDS"] = (
        ("nDS"),
        [
            (ds_era5["origin_datetime"] - np.datetime64("1970-01-01T00:00"))
            / np.timedelta64(1, "s")
        ],
        kpt_variables["timDS"],
    )
    ds_kpt["latDS"] = (("nDS"), [ds_era5["origin_lat"]], kpt_variables["latDS"])
    ds_kpt["lonDS"] = (("nDS"), [ds_era5["origin_lon"]], kpt_variables["lonDS"])
    # Change order of data, to confirm to other kpt input
    ds_kpt = ds_kpt.sortby("nlev", ascending=True)
    ds_kpt = ds_kpt.sortby("nlevp1", ascending=True)
    for var in kpt_variables:
        if var not in ds_kpt:
            print(var + " is missing in the kpt formatted output")
    # Needs improvement still
    kpt_dict = {
        "campaign": metadata.campaign,
        "flight": metadata.case,
        "date": ds_era5["origin_datetime"].values.astype("str"),
        "source": "ERA5",
        "source_domain": "NEEDS ADDING",
        "source_grid": "grid0.1x0.1",
        "source_latsamp": ds_era5.sampling_method[1],
        "source_lonsamp": ds_era5.sampling_method[1],
        "creator": metadata.author + " with https://github.com/EUREC4A-UK/lagtraj",
        "created": datetime.datetime.now().isoformat(),
        "wilting_point": 0.1715,
        "field_capacity": 0.32275,
        "t_skin_correct": "Skin temperature has been corrected "
        "by 1.000000. Motivation: value from IFS is actually "
        "the open SST, which is lower than the skin temperature.",
    }
    ds_kpt.attrs.update(**kpt_dict)
    return ds_kpt


def dephy_from_era5(ds_era5, da_levels, parameters, metadata):
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
    rh = rh_dephy(ds_dephy["temp"], ds_dephy["pressure"], ds_dephy["qt"])
    ds_dephy["rh"] = init_field_dephy(rh.values[:, :, 0, 0], "rh")

    def nudging_inv_time_prof(nudging_parameters, variable):
        height_array = ds_dephy["height_forc"].values
        nudging_method = nudging_parameters.method
        nudging_time = nudging_parameters.time
        nudging_height = nudging_parameters.height
        nudging_transition = nudging_parameters.transition
        if nudging_method == None:
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


def export(file_path, ds_conversion, nc_format=None):
    Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    encoding = validation.build_valid_encoding(ds=ds_conversion)
    ds_conversion.to_netcdf(file_path, encoding=encoding, format=nc_format)


# kpt variable : era5 variable
# (we loop over kpt variables here)
kpt_from_era5_variables = {
    "zf": "height_h_local",
    "zh": "height_h_local",
    "ps": "sp_mean",
    "pres": "p_h_mean",
    "presh": "p_h_mean",
    "u": "u_mean",
    "v": "v_mean",
    "t": "t_mean",
    "q": "q_mean",
    "ql": "clwc_mean",
    "qi": "ciwc_mean",
    "cloud_fraction": "cc_mean",
    "omega": "w_pressure_corr_mean",
    "o3": "o3_mean",
    "t_local": "t_local",
    "q_local": "q_local",
    "ql_local": "clwc_local",
    "qi_local": "ciwc_local",
    "u_local": "u_local",
    "v_local": "v_local",
    "cc_local": "cc_local",
    "tadv": "dtdt_adv",
    "qadv": "dqdt_adv",
    "uadv": "dudt_adv",
    "vadv": "dvdt_adv",
    "ug": "u_g",
    "vg": "v_g",
    "tladv": "dt_ldt_adv",
    "qladv": "dclwcdt_adv",
    "qiadv": "dciwcdt_adv",
    "ccadv": "dccdt_adv",
    "lat_grid": "lat",
    "lon_grid": "lon",
    "u_traj": "u_traj",
    "v_traj": "v_traj",
    "albedo": "fal_mean",
    "mom_rough": "fsr_mean",
    "t_snow": "tsn_mean",
    "albedo_snow": "asn_mean",
    "density_snow": "rsn_mean",
    "t_sea_ice": "istl1_mean",
    "open_sst": "sst_mean",
    "high_veg_cover": "cvh_mean",
    "low_veg_cover": "cvl_mean",
    "sea_ice_frct": "siconc_mean",
    "sdor": "sdor_mean",
    "isor": "isor_mean",
    "anor": "anor_mean",
}

# era5 units : kpt units
# we replace era5 units here
era5_to_kpt_units = {
    "m s**-1": "m/s",
    "1": "0-1",
    "degrees_north": "degrees North",
    "degrees_east": "degrees East",
    "metres": "m",
    "kg kg**-1": "kg/kg",
    "Pa s**-1": "Pa/s",
    "K s**-1": "K/s",
    "kg kg**-1 s**-1": "kg/kg/s",
    "m s**-1 s**-1": "m/s2",
    "kg m**-3": "kg/m3",
    "s**-1": "1/s",
    "W m**-2": "W/m2",
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


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("forcing")
    argparser.add_argument("conversion")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument("--debug", default=False, action="store_true")
    args = argparser.parse_args()

    try:
        ds_forcing = load_forcing_data(
            root_data_path=args.data_path, forcing_name=args.forcing
        )
    except FileNotFoundError:
        raise Exception(
            f"The output file for forcing `{args.forcing}`"
            " couldn't be found. Please create the trajectory by running: \n"
            f"    python -m lagtraj.forcings.create {args.forcing}\n"
            "and then run the forcing creation again"
        )

    conversion_defn = load.load_definition(
        root_data_path=args.data_path, conversion_name=args.conversion
    )

    if conversion_defn.levels.method is None or conversion_defn.levels.method == "copy":
        da_levels = ds_forcing["level"]
    else:
        da_levels = make_levels(
            method=conversion_defn.levels.method,
            n_levels=conversion_defn.levels.n_levels,
            z_top=conversion_defn.levels.z_top,
            dz_min=conversion_defn.levels.dz_min,
        )

    with optional_debugging(args.debug):
        if conversion_defn.export_format == "kpt":
            ds_conversion = kpt_from_era5(
                ds_forcing,
                da_levels,
                parameters=conversion_defn.parameters,
                metadata=conversion_defn.metadata,
            )
            nc_format = None
        elif conversion_defn.export_format == "dephy":
            ds_conversion = dephy_from_era5(
                ds_forcing,
                da_levels,
                parameters=conversion_defn.parameters,
                metadata=conversion_defn.metadata,
            )
            nc_format = "NETCDF3_CLASSIC"
        else:
            raise NotImplementedError(format)

    output_file_path = build_conversion_data_path(
        root_data_path=args.data_path,
        forcing_name=args.forcing,
        conversion_name=conversion_defn.name,
    )

    export(ds_conversion=ds_conversion, file_path=output_file_path, nc_format=nc_format)

    print("Wrote forcing file to `{}`".format(output_file_path))


if __name__ == "__main__":
    main()
