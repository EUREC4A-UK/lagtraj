"""This module contains:
    1. a list of the KPT variables as a dictionary
    2. the mapping between the internal lagtraj variables (era5) and KPT
    3. and a function to process the conversion
"""
import datetime

import numpy as np
import xarray as xr

from ....domain.sources.era5.constants import rg
from ....utils.interpolation.methods import (
    central_estimate,
    cos_transition,
    steffen_1d_no_ep_time,
)

kpt_attributes = {
    "lat": {"units": "degrees North", "long_name": "latitude"},
    "lon": {"units": "degrees East", "long_name": "longitude"},
    "zf": {"units": "m", "long_name": "full level height"},
    "zh": {"units": "m", "long_name": "half level height"},
    "ps": {"units": "Pa", "long_name": "surface pressure"},
    "pres": {"units": "Pa", "long_name": "full level pressure"},
    "presh": {"units": "Pa", "long_name": "half level pressure"},
    "u": {"units": "m/s", "long_name": "zonal wind (domain averaged)"},
    "v": {"units": "m/s", "long_name": "meridional wind (domain averaged)"},
    "t": {"units": "K", "long_name": "temperature (domain averaged)"},
    "q": {"units": "kg/kg", "long_name": "water vapor mixing ratio (domain averaged)"},
    "ql": {
        "units": "kg/kg",
        "long_name": "liquid water mixing ratio (domain averaged)",
    },
    "qi": {"units": "kg/kg", "long_name": "ice water mixing ratio (domain averaged)"},
    "cloud_fraction": {"units": "0-1", "long_name": "cloud fraction (domain averaged)"},
    "omega": {
        "units": "Pa/s",
        "long_name": "large-scale pressure velocity (domain averaged)",
    },
    "o3": {"units": "kg/kg", "long_name": "ozone mass mixing ratio (domain averaged)"},
    "t_local": {"units": "K", "long_name": "temperature (at domain midpoint)"},
    "q_local": {
        "units": "kg/kg",
        "long_name": "water vapor specific humidity (at domain midpoint)",
    },
    "ql_local": {
        "units": "kg/kg",
        "long_name": "liquid water specific humidity (at domain midpoint)",
    },
    "qi_local": {
        "units": "kg/kg",
        "long_name": "ice water specific humidity (at domain midpoint)",
    },
    "u_local": {"units": "m/s", "long_name": "zonal wind (at domain midpoint)"},
    "v_local": {"units": "m/s", "long_name": "meridional wind (at domain midpoint)"},
    "cc_local": {"units": "0-1", "long_name": "cloud fraction (at domain midpoint)"},
    "tadv": {
        "units": "K/s",
        "long_name": "tendency in temperature due to "
        "large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using velocity relative to wind "
        "on trajectory (u_traj,v_traj)",
    },
    "qadv": {
        "units": "kg/kg/s",
        "long_name": "tendency in water vapor due to large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using velocity relative to wind on trajectory (u_traj,v_traj)",
    },
    "uadv": {
        "units": "m/s2",
        "long_name": "tendency in zonal wind due to large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using velocity relative to wind on trajectory (u_traj,v_traj)",
    },
    "vadv": {
        "units": "m/s2",
        "long_name": "tendency in meridional wind due to large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using velocity relative to wind on trajectory (u_traj,v_traj)",
    },
    "ug": {
        "units": "m/s",
        "long_name": "geostrophic wind - zonal component",
        "info": "derived at pressure levels",
        "interpolation": "above 5 hPa the geostrophic wind is equal to the real wind",
    },
    "vg": {
        "units": "m/s",
        "long_name": "geostrophic wind -meridional component",
        "info": "derived at pressure levels",
        "interpolation": "above 5 hPa the geostrophic wind is equal to the real wind",
    },
    "tladv": {
        "units": "K/s",
        "long_name": "tendency in T_l due to large-scale horizontal advection",
    },
    "qladv": {
        "units": "kg/kg/s",
        "long_name": "tendency in liquid water spec hum due to large-scale horizontal advection",
    },
    "qiadv": {
        "units": "kg/kg/s",
        "long_name": "tendency in frozen water due to large-scale horizontal advection",
    },
    "ccadv": {
        "units": "1/s",
        "long_name": "tendency in cloud fraction due to large-scale horizontal advection",
    },
    "time_traj": {
        "units": "seconds since 1-1-1970 00:00",
        "long_name": "time values at trajectory waypoints",
    },
    "lat_traj": {
        "units": "degrees North",
        "long_name": "latitude of trajectory waypoints",
    },
    "lon_traj": {
        "units": "degrees East",
        "long_name": "longitude of trajectory waypoints",
    },
    "u_traj": {"units": "m/s", "long_name": "zonal wind at trajectory waypoints"},
    "v_traj": {"units": "m/s", "long_name": "meridional wind at trajectory waypoints"},
    "albedo": {"units": "0-1", "long_name": "albedo"},
    "mom_rough": {"units": "m", "long_name": "roughness length for momentum"},
    "heat_rough": {"units": "m", "long_name": "roughness length for heat"},
    "t_skin": {
        "units": "K",
        "long_name": "skin temperature",
    },
    "q_skin": {"units": "m of water", "long_name": "skin reservoir content"},
    "snow": {"units": "m, liquid equivalent", "long_name": "snow depth"},
    "t_snow": {"units": "K", "long_name": "snow temperature"},
    "albedo_snow": {"units": "0-1", "long_name": "snow albedo"},
    "density_snow": {"units": "kg/m3", "long_name": "snow density"},
    "sfc_sens_flx": {"units": "W/m2", "long_name": "surface sensible heat flux"},
    "sfc_lat_flx": {"units": "W/m2", "long_name": "surface latent heat flux"},
    "h_soil": {"units": "m", "long_name": "soil layer thickness"},
    "t_soil": {"units": "K", "long_name": "soil layer temperature"},
    "q_soil": {"units": "m3/m3", "long_name": "soil moisture"},
    "lsm": {"units": "-", "long_name": "land sea mask"},
    "t_sea_ice": {"units": "K", "long_name": "sea ice temperature"},
    "open_sst": {"units": "K", "long_name": "open sea surface temperature"},
    "orog": {"units": "m", "long_name": "orography - surface height"},
    "DS": {"long_name": "label of trajectory reference point"},
    "timDS": {
        "units": "seconds since 1-1-1970 00:00",
        "long_name": "time at trajectory reference point",
        "info": "the reference point is the space-time "
        "coordinate from which the trajectory is calculated",
    },
    "latDS": {
        "units": "degrees North",
        "long_name": "latitude at trajectory reference point",
        "info": "the reference point is the space-time coordinate "
        "from which the trajectory is calculated",
    },
    "lonDS": {
        "units": "degrees East",
        "long_name": "longitude at trajectory reference point",
        "info": "the reference point is the space-time coordinate "
        "from which the trajectory is calculated",
    },
    "lat_grid": {
        "units": "degrees North",
        "long_name": "latitude of closest IFS gridpoint",
    },
    "lon_grid": {
        "units": "degrees East",
        "long_name": "longitude of closest IFS gridpoint",
    },
    "p_traj": {
        "units": "hPa",
        "long_name": "pressure level at which trajectory was calculated",
    },
    "sv": {"units": "whatever", "long_name": "tracers"},
    "fradSWnet": {"units": "W/m2", "long_name": "radiative flux - net short wave"},
    "fradLWnet": {"units": "W/m2", "long_name": "radiative flux - net long wave"},
    "high_veg_type": {"units": "-", "long_name": "high vegetation type"},
    "low_veg_type": {"units": "-", "long_name": "low vegetation type"},
    "high_veg_cover": {"units": "0-1", "long_name": "high vegetation cover"},
    "low_veg_cover": {"units": "0-1", "long_name": "low vegetation cover"},
    "high_veg_lai": {"units": "-", "long_name": "leaf area index of high vegetation"},
    "low_veg_lai": {"units": "-", "long_name": "leaf area index of low vegetation"},
    "sea_ice_frct": {"units": "0-1", "long_name": "sea ice fraction"},
    "sdor": {
        "units": "m",
        "long_name": "subgrid-scale orography - standard deviation",
    },
    "isor": {"units": "0-1", "long_name": "subgrid-scale orography - anisotropy"},
    "anor": {
        "units": "radians",
        "long_name": "subgrid-scale orography - orientation/angle of steepest gradient",
    },
    "slor": {"units": "m/m", "long_name": "subgrid-scale orography - mean slope"},
    "msnswrf": {
        "units": "W/m2",
        "long_name": "Mean surface net short-wave radiation flux",
    },
    "msnlwrf": {
        "units": "W/m2",
        "long_name": "Mean surface net long-wave radiation flux",
    },
    "mtnswrf": {
        "units": "W/m2",
        "long_name": "Mean top net short-wave radiation flux",
    },
    "mtnlwrf": {"units": "W/m2", "long_name": "Mean top net long-wave radiation flux"},
    "mtnswrfcs": {
        "units": "W/m2",
        "long_name": "Mean top net short-wave radiation flux, clear sky",
    },
    "mtnlwrfcs": {
        "units": "W/m2",
        "long_name": "Mean top net long-wave radiation flux, clear sky",
    },
    "msnswrfcs": {
        "units": "W/m2",
        "long_name": "Mean surface net short-wave radiation flux, clear sky",
    },
    "msnlwrfcs": {
        "units": "W/m2",
        "long_name": "Mean surface net long-wave radiation flux, clear sky",
    },
    "mtdwswrf": {
        "units": "W/m2",
        "long_name": "Mean top downward short-wave radiation flux",
    },
    "second": {"units": "s", "long_name": "second of day"},
    "date": {"units": "days", "long_name": "date in yyyymmdd format"},
}


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
    "lat": "lat",
    "lon": "lon",
    "lat_traj": "lat",
    "lon_traj": "lon",
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

# create a special mapping to allow for variations in units of data coming from
# ERA5, since changes have happened over time
ALLOWED_UNIT_VARIATIONS = dict(
    # `sdor` (standard deviation of orography) previously had units of `0-1`
    # but this was corrected to `m` in ecCodes as of version v2.22.0 (May
    # 2021). We allow for the old units and always set to the correct value of
    # `m` for kpt output (see https://github.com/EUREC4A-UK/lagtraj/pull/169)
    sdor=dict(valid_era5_units=("0-1", "m"), kpt_units="m")
)


def from_era5(ds_era5, da_levels, parameters, metadata):
    """Obtain a kpt input file from era5 variable set at high resolution"""
    # Put full levels midway between half-levels, I think this is consistent with DALES
    # Reverse order of data, to confirm to other kpt input
    kpt_half_level_array = da_levels.values
    kpt_full_level_array = 0.5 * (kpt_half_level_array[:-1] + kpt_half_level_array[1:])
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
    nDS_coord = {
        "nDS": (
            "nDS",
            [0],
            {},
        )
    }
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
        if not unit_guess == kpt_attributes[variable]["units"]:
            kpt_units = kpt_attributes[variable]["units"]
            if variable in ALLOWED_UNIT_VARIATIONS:
                unit_variation = ALLOWED_UNIT_VARIATIONS[variable]
                era5_has_valid_units = unit_guess in unit_variation["valid_era5_units"]
                kpt_units_valid = kpt_units == unit_variation["kpt_units"]
                if era5_has_valid_units and kpt_units_valid:
                    # the units we want to use for KPT and the ones that we say
                    # that ERA5 is actually in
                    pass
                elif not era5_has_valid_units:
                    raise Exception(
                        f"The variable `{variable}` has changed units in the ERA5 source files"
                        " over time, but the ERA5 units don't match the expected value values"
                        f" {', '.join(unit_variation['valid_era5_units'])}"
                    )
                else:
                    raise Exception(
                        f"The variable `{variable}` has changed units in the ERA5 source files"
                        f" over time, but the units requested for KPT output `{kpt_units}`"
                        f" does not match that selected for the ERA5 variable {kpt_units_valid}"
                    )
            else:
                raise Exception(
                    f"Incompatible units between ERA5 and kpt for variable `{variable}`. "
                    "Please fix using the era5_to_kpt_units dictionary: "
                    f"ERA converted variable is {unit_guess}, kpt variable is {kpt_units}"
                )
        # single level variable
        if np.ndim(da_era5.values) == 1:
            ds_kpt[variable] = da_era5
            ds_kpt[variable].attrs.update(kpt_attributes[variable])
        # half level variable
        elif variable in ["zh", "presh"]:
            da_era5_on_half_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, kpt_half_level_array
            )
            ds_kpt[variable] = (
                ("time", "nlevp1"),
                da_era5_on_half_levels,
                kpt_attributes[variable],
            )
        # full level variable
        else:
            da_era5_on_full_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, kpt_full_level_array
            )
            ds_kpt[variable] = (
                ("time", "nlev"),
                da_era5_on_full_levels,
                kpt_attributes[variable],
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
        ds_kpt[variable].attrs.update(**kpt_attributes[variable])
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
            kpt_attributes[variable],
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
        kpt_attributes["q_soil"],
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
        kpt_attributes["t_soil"],
    )
    # Soil thickness: combine levels
    h_soil = np.array([0.07, 0.21, 0.72, 1.89])
    ds_kpt["h_soil"] = (
        ("nlevs"),
        h_soil,
        kpt_attributes["h_soil"],
    )
    # Orography: derive from surface geopotential
    ds_kpt["orog"] = ds_era5["z_mean"] / rg
    ds_kpt["orog"].attrs.update(**kpt_attributes["orog"])
    # Heat roughness, derive from "flsr" variable
    ds_kpt["heat_rough"] = np.exp(ds_era5["flsr_mean"])
    ds_kpt["heat_rough"].attrs.update(**kpt_attributes["heat_rough"])
    ds_kpt["t_skin"] = ds_era5["skt_mean"]
    ds_kpt["t_skin"].attrs.update(**kpt_attributes["t_skin"])
    # Surface fluxes: obtain from time mean in ERA data, do not change sign
    sfc_sens_flx = central_estimate(ds_era5["msshf_mean"].values)
    ds_kpt["sfc_sens_flx"] = (
        ("time"),
        sfc_sens_flx,
        kpt_attributes["sfc_sens_flx"],
    )
    sfc_lat_flx = central_estimate(ds_era5["mslhf_mean"].values)
    ds_kpt["sfc_lat_flx"] = (
        ("time"),
        sfc_lat_flx,
        kpt_attributes["sfc_lat_flx"],
    )
    # Final checks: are all variables present?
    ds_kpt["time_traj"] = (
        ds_era5["time"] - np.datetime64("1970-01-01T00:00")
    ) / np.timedelta64(1, "s")
    ds_kpt["time_traj"].attrs.update(**kpt_attributes["time_traj"])
    ds_kpt["DS"] = (("nDS"), ["Trajectory origin"], kpt_attributes["DS"])
    ds_kpt["timDS"] = (
        ("nDS"),
        [
            (ds_era5["origin_datetime"] - np.datetime64("1970-01-01T00:00"))
            / np.timedelta64(1, "s")
        ],
        kpt_attributes["timDS"],
    )
    ds_kpt["latDS"] = (("nDS"), [ds_era5["origin_lat"]], kpt_attributes["latDS"])
    ds_kpt["lonDS"] = (("nDS"), [ds_era5["origin_lon"]], kpt_attributes["lonDS"])
    # Change order of data, to confirm to other kpt input
    ds_kpt = ds_kpt.sortby("nlev", ascending=True)
    ds_kpt = ds_kpt.sortby("nlevp1", ascending=True)
    ds_kpt["second"] = ds_kpt["time_traj"] % 86400
    ds_kpt["date"] = (
        ("time"),
        (ds_kpt["time"].dt.strftime("%Y%m%d")).values.astype("int"),
        kpt_attributes["date"],
    )
    for var in kpt_attributes:
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
        # TODO: (Leif) these need adding, where should they come from?
        # "source_latsamp": ds_era5.sampling_method[1],
        # "source_lonsamp": ds_era5.sampling_method[1],
        "creator": f"{metadata.author} with https://github.com/EUREC4A-UK/lagtraj",
        "created": datetime.datetime.now().isoformat(),
        "wilting_point": 0.1715,
        "field_capacity": 0.32275,
    }
    ds_kpt.attrs.update(**kpt_dict)
    if parameters.inversion_nudging is not None:
        if parameters.inversion_nudging in [0, 1]:
            ds_inversion = {
                "inversion_nudging": parameters.inversion_nudging,
                "inversion_nudging_height_above": parameters.inversion_nudging_height_above,
                "inversion_nudging_transition": parameters.inversion_nudging_transition,
                "inversion_nudging_time": parameters.inversion_nudging_time,
            }
        else:
            raise NotImplementedError(
                f"Inversion nudging option `{parameters.inversion_nudging}` not implemented"
            )
        ds_kpt.attrs.update(**ds_inversion)
    # Correct geostropic winds and wind tendencies at high levels
    if parameters.wind_at_high_levels_correction is None:
        # Use sensible default values
        wind_at_high_levels_correction = 1
        wind_at_high_levels_correction_pressure_above = 500.0  # Pa, not hPa!
        wind_at_high_levels_correction_transition = 200.0
    elif parameters.wind_at_high_levels_correction in [0, 1]:
        wind_at_high_levels_correction = parameters.wind_at_high_levels_correction
        wind_at_high_levels_correction_pressure_above = (
            parameters.wind_at_high_levels_correction_pressure_above
        )  # Pa, not hPa!
        wind_at_high_levels_correction_transition = (
            parameters.wind_at_high_levels_correction_transition
        )
    else:
        raise NotImplementedError(
            f"Wind at high level correction option `{parameters.wind_at_high_levels_correction}` not implemented"
        )
    ds_wind_at_high_levels = {
        "wind_at_high_levels_correction": wind_at_high_levels_correction,
        "wind_at_high_levels_correction_pressure_above": wind_at_high_levels_correction_pressure_above,
        "wind_at_high_levels_correction_transition": wind_at_high_levels_correction_transition,
    }
    ds_kpt.attrs.update(**ds_wind_at_high_levels)
    pressure_array = ds_kpt["pres"].values
    wind_at_high_levels_correction_factor = cos_transition(
        pressure_array,
        wind_at_high_levels_correction_pressure_above
        + 0.5 * wind_at_high_levels_correction_transition,
        wind_at_high_levels_correction_pressure_above
        - 0.5 * wind_at_high_levels_correction_transition,
    )
    # Set ug and vg equal to actual (nudging) wind at very high levels
    # Remove advection tendencies at high levels
    if wind_at_high_levels_correction == 1:
        ds_kpt["ug"] = ds_kpt["ug"] * wind_at_high_levels_correction_factor[:] + ds_kpt[
            "u"
        ] * (1.0 - wind_at_high_levels_correction_factor)
        ds_kpt["vg"] = ds_kpt["vg"] * wind_at_high_levels_correction_factor + ds_kpt[
            "v"
        ] * (1.0 - wind_at_high_levels_correction_factor)
        ds_kpt["uadv"] = ds_kpt["uadv"] * wind_at_high_levels_correction_factor
        ds_kpt["vadv"] = ds_kpt["vadv"] * wind_at_high_levels_correction_factor
        ds_kpt["tadv"] = ds_kpt["tadv"] * wind_at_high_levels_correction_factor
        ds_kpt["qadv"] = ds_kpt["qadv"] * wind_at_high_levels_correction_factor
        ds_kpt["tladv"] = ds_kpt["tladv"] * wind_at_high_levels_correction_factor
        ds_kpt["qladv"] = ds_kpt["qladv"] * wind_at_high_levels_correction_factor
        ds_kpt["qiadv"] = ds_kpt["qiadv"] * wind_at_high_levels_correction_factor
        ds_kpt["ccadv"] = ds_kpt["ccadv"] * wind_at_high_levels_correction_factor
    return ds_kpt
