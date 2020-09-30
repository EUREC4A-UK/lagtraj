import xarray as xr
from tqdm import tqdm
from pathlib import Path
import datetime
import numpy as np

from .. import DEFAULT_ROOT_DATA_PATH
from ..forcings.load import load_data as load_forcing_data
from . import load, build_conversion_data_path
from ..utils import optional_debugging, validation
from ..utils.xarray import create_attributes_dictionary  # Needs using?
from ..utils.thermo import rh_hightune
from ..utils.interpolation import steffen_1d_no_ep_time, central_estimate
from .racmo import racmo_variables
from .hightune import hightune_variables
from ..domain.sources.era5.constants import rg, cp, rlv

# TODO:
# - Only interpolate if needed?
# - Check forcing has the right format?
# Fix time formats
# Fix missing variables


def racmo_from_era5(ds_era5, da_levels):
    """Obtain a racmo input file from era5 variable set at high resolution"""
    # Put full levels midway between half-levels, I think this is consistent with DALES
    # Reverse order of data, to confirm to other RACMO input
    racmo_half_level_array = da_levels.values
    racmo_full_level_array = 0.5 * (
        racmo_half_level_array[:-1] + racmo_half_level_array[1:]
    )
    racmo_half_level_coord = {
        "nlevp1": (
            "nlevp1",
            (np.arange(len(racmo_half_level_array)) + 1.0)[::-1],
            {"long_name": "model half levels"},
        )
    }
    racmo_full_level_coord = {
        "nlev": (
            "nlev",
            (np.arange(len(racmo_full_level_array)) + 1.0)[::-1],
            {"long_name": "model full levels"},
        )
    }
    nDS_coord = {"nDS": ("nDS", [0], {},)}
    racmo_soil_coord = {
        "nlevs": ("nlevs", np.arange(4) + 1.0, {"long_name": "soil levels"})
    }
    ds_racmo = xr.Dataset(
        coords={
            "time": ds_era5.time,
            **racmo_full_level_coord,
            **racmo_half_level_coord,
            **racmo_soil_coord,
            **nDS_coord,
        }
    )
    # Variables from dictionary
    # Including unit checks
    for variable in racmo_from_era5_variables:
        era5_var = racmo_from_era5_variables[variable]
        da_era5 = ds_era5[era5_var]
        # perform units check
        unit_guess = era5_to_racmo_units.get(da_era5.units, da_era5.units)
        if not (unit_guess == racmo_variables[variable]["units"]):
            except_str = (
                "Incompatible units between ERA5 and RACMO for variable "
                + variable
                + ". Please fix using the fix_era5_to_racmo_units dictionary: ERA converted variable is "
                + unit_guess
                + ", RACMO variable is "
                + racmo_variables[variable]["units"]
            )
            raise Exception(except_str)
        # single level variable
        if np.ndim(da_era5.values) == 1:
            ds_racmo[variable] = (("time"), da_era5, racmo_variables[variable])
        # half level variable
        elif variable in ["zh", "presh"]:
            da_era5_on_half_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, racmo_half_level_array
            )
            ds_racmo[variable] = (
                ("time", "nlevp1"),
                da_era5_on_half_levels,
                racmo_variables[variable],
            )
        # full level variable
        else:
            da_era5_on_full_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, racmo_full_level_array
            )
            ds_racmo[variable] = (
                ("time", "nlev"),
                da_era5_on_full_levels,
                racmo_variables[variable],
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
        ds_racmo[variable] = ds_era5[variables_to_manually_add[variable]]
        ds_racmo[variable].attrs.update(**racmo_variables[variable])
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
        ds_racmo[variable] = (
            ("time"),
            this_central_estimate,
            racmo_variables[variable],
        )
    # Soil moisture: combine levels
    swvl1 = ds_era5["swvl1_mean"].values
    swvl2 = ds_era5["swvl2_mean"].values
    swvl3 = ds_era5["swvl3_mean"].values
    swvl4 = ds_era5["swvl4_mean"].values
    q_soil = np.stack((swvl1, swvl2, swvl3, swvl4), axis=-1)
    ds_racmo["q_soil"] = (
        ("time", "nlevs"),
        q_soil,
        racmo_variables["q_soil"],
    )
    # Soil temperature: combine levels
    stl1 = ds_era5["stl1_mean"].values
    stl2 = ds_era5["stl2_mean"].values
    stl3 = ds_era5["stl3_mean"].values
    stl4 = ds_era5["stl4_mean"].values
    t_soil = np.stack((stl1, stl2, stl3, stl4), axis=-1)
    ds_racmo["t_soil"] = (
        ("time", "nlevs"),
        t_soil,
        racmo_variables["t_soil"],
    )
    # Soil thickness: combine levels
    h_soil = np.array([0.07, 0.21, 0.72, 1.89])
    ds_racmo["h_soil"] = (
        ("nlevs"),
        h_soil,
        racmo_variables["h_soil"],
    )
    # Orography: derive from surface geopotential
    ds_racmo["orog"] = ds_era5["z_mean"] / rg
    ds_racmo["orog"].attrs.update(**racmo_variables["orog"])
    # Heat roughness, derive from "flsr" variable
    ds_racmo["heat_rough"] = np.exp(ds_era5["flsr_mean"])
    ds_racmo["heat_rough"].attrs.update(**racmo_variables["heat_rough"])
    # Apply correction to t_skin (see output files)
    ds_racmo["t_skin"] = ds_era5["stl1_mean"] + 1.0
    ds_racmo["t_skin"].attrs.update(**racmo_variables["t_skin"])
    # Surface fluxes: obtain from time mean in ERA data, do not change sign
    sfc_sens_flx = central_estimate(ds_era5["msshf_mean"].values)
    ds_racmo["sfc_sens_flx"] = (
        ("time"),
        sfc_sens_flx,
        racmo_variables["sfc_sens_flx"],
    )
    sfc_lat_flx = central_estimate(ds_era5["mslhf_mean"].values)
    ds_racmo["sfc_lat_flx"] = (
        ("time"),
        sfc_lat_flx,
        racmo_variables["sfc_lat_flx"],
    )
    # Final checks: are all variables present?
    ds_racmo["time_traj"] = (
        ds_era5["time"] - np.datetime64("1970-01-01T00:00")
    ) / np.timedelta64(1, "s")
    ds_racmo["time_traj"].attrs.update(**racmo_variables["time_traj"])
    ds_racmo["DS"] = (("nDS"), ["Trajectory origin"], racmo_variables["DS"])
    ds_racmo["timDS"] = (
        ("nDS"),
        [
            (ds_era5["origin_datetime"] - np.datetime64("1970-01-01T00:00"))
            / np.timedelta64(1, "s")
        ],
        racmo_variables["timDS"],
    )
    ds_racmo["latDS"] = (("nDS"), [ds_era5["origin_lat"]], racmo_variables["latDS"])
    ds_racmo["lonDS"] = (("nDS"), [ds_era5["origin_lon"]], racmo_variables["lonDS"])
    # Change order of data, to confirm to other RACMO input
    ds_racmo = ds_racmo.sortby("nlev", ascending=True)
    ds_racmo = ds_racmo.sortby("nlevp1", ascending=True)
    for var in racmo_variables:
        if var not in ds_racmo:
            print(var + " is missing in the RACMO formatted output")
    # Needs improvement still
    racmo_dict = {
        "campaign": "NEEDS ADDING",
        "flight": "NEEDS ADDING",
        "date": ds_era5["origin_datetime"].values.astype("str"),
        "source": "ERA5",
        "source_domain": "NEEDS ADDING",
        "source_grid": "grid0.1x0.1",
        "source_latsamp": ds_era5.sampling_method[1],
        "source_lonsamp": ds_era5.sampling_method[1],
        "creator": "https://github.com/EUREC4A-UK/lagtraj",
        "created": datetime.datetime.now().isoformat(),
        "wilting_point": 0.1715,
        "field_capacity": 0.32275,
        "t_skin_correct": "Skin temperature has been corrected by 1.000000. Motivation: value from IFS is actually the open SST, which is lower than the skin temperature.",
    }
    ds_racmo.attrs.update(**racmo_dict)
    return ds_racmo


# VARIABLES CURRENTLY LEFT OUT
# 't0' : {'long_name': 'Initial time'}
def hightune_from_era5(ds_era5, da_levels):
    """Obtain a hightune input file from era5 variable set at high resolution"""
    hightune_level_array = da_levels.values
    hightune_level_coord = {
        "lev": ("lev", da_levels.values, {"long_name": "altitude", "units": "m"},)
    }
    ds_hightune = xr.Dataset(coords={"time": ds_era5.time, **hightune_level_coord,})
    # Variables from dictionary
    # Including unit checks
    for variable in hightune_from_era5_variables:
        era5_var = hightune_from_era5_variables[variable]
        da_era5 = ds_era5[era5_var]
        # perform units check
        unit_guess = era5_to_hightune_units.get(da_era5.units, da_era5.units)
        if not (unit_guess == hightune_variables[variable]["units"]):
            except_str = (
                "Incompatible units between ERA5 and hightune for variable "
                + variable
                + ". Please fix using the fix_era5_to_hightune_units dictionary: ERA converted variable is "
                + unit_guess
                + ", hightune variable is "
                + hightune_variables[variable]["units"]
            )
            raise Exception(except_str)
        if np.ndim(da_era5.values) == 1:
            ds_hightune[variable] = (("time"), da_era5, hightune_variables[variable])
        else:
            da_era5_on_hightune_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, hightune_level_array
            )
            ds_hightune[variable] = (
                ("time", "lev"),
                da_era5_on_hightune_levels,
                hightune_variables[variable],
            )
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
        ds_hightune[variable] = (
            ("time"),
            this_central_estimate,
            hightune_variables[variable],
        )
    # TKE, set to zero
    ds_hightune["tke"] = 0.0 * (ds_era5["u_mean"] * ds_era5["u_mean"])
    ds_hightune["tke"].attrs.update(**hightune_variables["tke"])
    # Heat roughness, derive from "flsr" variable
    ds_hightune["z0h"] = np.exp(ds_era5["flsr_mean"])
    ds_hightune["z0h"].attrs.update(**hightune_variables["z0h"])
    # Include same t_skin correction used for DALES, may need further work
    ds_hightune["ts"] = ds_era5["stl1_mean"] + 1.0
    ds_hightune["ts"].attrs.update(**hightune_variables["ts"])
    # Surface fluxes: obtain from time mean in ERA data, change sign for hightune!
    sfc_sens_flx = -central_estimate(ds_era5["msshf_mean"].values)
    ds_hightune["sfc_sens_flx"] = (
        ("time"),
        sfc_sens_flx,
        hightune_variables["sfc_sens_flx"],
    )
    sfc_lat_flx = -central_estimate(ds_era5["mslhf_mean"].values)
    ds_hightune["sfc_lat_flx"] = (
        ("time"),
        sfc_lat_flx,
        hightune_variables["sfc_lat_flx"],
    )
    wpthetap = (sfc_sens_flx / (cp * ds_era5["rho_mean"].sel(level=0.0))) * (
        ds_era5["theta_mean"].sel(level=0.0) / ds_era5["t_mean"].sel(level=0.0)
    )
    ds_hightune["wpthetap"] = (
        ("time"),
        wpthetap,
        hightune_variables["wpthetap"],
    )
    wpqvp = sfc_lat_flx / (rlv * ds_era5["rho_mean"].sel(level=0.0))
    ds_hightune["wpqvp"] = (
        ("time"),
        wpqvp,
        hightune_variables["wpqvp"],
    )
    wpqtp = wpqvp
    ds_hightune["wpqtp"] = (
        ("time"),
        wpqtp,
        hightune_variables["wpqtp"],
    )
    # Ratio of fluxes (mixing ratio vs. specific humidity) is same for all fluxes
    moisture_ratio = ds_era5["r_t_local"].sel(level=0.0) / ds_era5["q_t_local"].sel(
        level=0.0
    )
    wprvp = wpqvp * moisture_ratio
    ds_hightune["wprvp"] = (
        ("time"),
        wprvp,
        hightune_variables["wprvp"],
    )
    wprtp = wpqtp * moisture_ratio
    ds_hightune["wprtp"] = (
        ("time"),
        sfc_lat_flx,
        hightune_variables["wprtp"],
    )
    ds_hightune["rh"] = rh_hightune(
        ds_hightune["temp"], ds_hightune["pressure"], ds_hightune["qt"],
    )
    ds_hightune["rh"].attrs.update(**hightune_variables["rh"])
    # Final checks: are all variables present?
    for var in hightune_variables:
        if var not in ds_hightune:
            print(var + " is missing in the hightune formatted output")
    # Needs improvement
    hightune_dictionary = {
        "Conventions": "CF-1.0",
        "comment": "Forcing and initial conditions for Lagrangian case",
        "reference": "NEEDS ADDING",
        "author": "NEEDS ADDING",
        "modifications": "NEEDS ADDING",
        "case": "NEEDS ADDING",
        "script": "https://github.com/EUREC4A-UK/lagtraj",
        "startDate": ds_hightune["time"][0].values.astype("str"),
        "endDate": ds_hightune["time"][-1].values.astype("str"),
        "tadv": 0,
        "tadvh": 1,
        "tadvv": 0,
        "rad_temp": 1,
        "qvadv": 0,
        "qvadvh": 1,
        "qvadvv": 0,
        "forc_omega": 0,
        "forc_w": 1,
        "forc_geo": 1,
        "nudging_u": 0,
        "nudging_v": 0,
        "nudging_t": 0,
        "nudging_q": 0,
        "zorog": 0.0,
        "surfaceType": "ocean",
        "surfaceForcing": "ts",
        "surfaceForcingWind": "z0_lagtraj",
    }
    ds_hightune.attrs.update(hightune_dictionary)
    return ds_hightune


def export(file_path, ds_conversion):
    Path(file_path).parent.mkdir(parents=True, exist_ok=True)
    encoding = validation.build_valid_encoding(ds=ds_conversion)
    ds_conversion.to_netcdf(file_path, encoding=encoding)


# racmo variable : era5 variable
# (we loop over racmo variables here)
racmo_from_era5_variables = {
    "lat": "lat",
    "lon": "lon",
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
    "lat_traj": "lat",
    "lon_traj": "lon",
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

# era5 units : racmo units
# we replace era5 units here
era5_to_racmo_units = {
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

# era5 units : hightune units
# we replace era5 units here
era5_to_hightune_units = {
    "m s**-1": "m s-1",
    "metres": "m",
    "kg kg**-1": "kg kg-1",
    "K s**-1": "K s-1",
    "kg kg**-1 s**-1": "kg kg-1 s-1",
    "W m**-2": "W m-2",
}


# hightune variable : era5 variable
# we loop over hightune variables here
hightune_from_era5_variables = {
    "lat": "lat",
    "lon": "lon",
    "ps": "sp_mean",
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
    "ps_forc": "sp_mean",
    "height_forc": "height_h_mean",
    "pressure_forc": "p_h_mean",
    "ug": "u_g",
    "vg": "v_g",
    "temp_adv": "dtdt_adv",
    "theta_adv": "dthetadt_adv",
    "thetal_adv": "dtheta_ldt_adv",
    "qv_adv": "dqdt_adv",
    "qt_adv": "dq_tdt_adv",
    "rv_adv": "dr_vdt_adv",
    "rt_adv": "dr_tdt_adv",
    "w": "w_corr_mean",
    "temp_nudging": "t_mean",
    "theta_nudging": "theta_mean",
    "thetal_nudging": "theta_l_mean",
    "qv_nudging": "q_mean",
    "qt_nudging": "q_t_mean",
    "rv_nudging": "r_v_mean",
    "rt_nudging": "r_t_mean",
    "u_nudging": "u_mean",
    "v_nudging": "v_mean",
    "z0m": "fsr_mean",
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
            f"The output file for forcing `{forcing_defn.name}`"
            " couldn't be found. Please create the trajectory by running: \n"
            f"    python -m lagtraj.forcings.create {forcing_defn.name}\n"
            "and then run the forcing creation again"
        )

    conversion_defn = load.load_definition(
        root_data_path=args.data_path, conversion_name=args.conversion
    )

    if conversion_defn.levels.method == None or conversion_defn.levels.method == "copy":
        da_levels = ds_forcing["level"]
    else:
        da_levels = make_levels(
            method=levels_definition.method,
            n_levels=levels_definition.n_levels,
            z_top=levels_definition.z_top,
            dz_min=levels_definition.dz_min,
        )

    with optional_debugging(args.debug):
        if conversion_defn.export_format == "racmo":
            ds_conversion = racmo_from_era5(ds_forcing, da_levels)
        elif conversion_defn.export_format == "hightune":
            ds_conversion = hightune_from_era5(ds_forcing, da_levels)
        else:
            raise NotImplementedError(format)

    output_file_path = build_conversion_data_path(
        root_data_path=args.data_path,
        forcing_name=args.forcing,
        conversion_name=conversion_defn.name,
    )

    export(ds_conversion=ds_conversion, file_path=output_file_path)

    print("Wrote forcing file to `{}`".format(output_file_path))


if __name__ == "__main__":
    main()
