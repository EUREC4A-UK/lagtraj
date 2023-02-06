"""This module contains:
    1. a list of the SAM variables as a dictionary
    2. the mapping between the internal lagtraj variables (era5) and SAM
    3. and a function to process the conversion
    4. a function to compute some time related variables for SAM.
    Notes: this module is created based on kpt.py because Peter Blossey's matlab conversion scripts
           take in kpt format output.
"""
import datetime

import numpy as np

# import pandas as pd
import scipy.interpolate as intp
import xarray as xr

# from ....domain.sources.era5.constants import rg
from ....utils.interpolation.methods import central_estimate, steffen_1d_no_ep_time
from ....utils.thermo import cpd, hrl, hrv

# following Peter Blossey's matlab script structure:
sam_attributes = {
    "Ps": {"units": "Pa", "long_name": "surface pressure"},
    "Ptend": {"units": "Pa/s", "long_name": "surface pressure tendency"},
    "Tg": {
        "units": "K",
        "long_name": "Surface (skin) Temperature (SST if over water)",
    },
    "shflx": {"units": "W/m2", "long_name": "surface sensible heat flux"},
    "lhflx": {"units": "W/m2", "long_name": "surface latent heat flux"},
    "Tsair": {
        "units": "K",
        "long_name": "Surface Air Temperature (extrapolated to z=0)",
    },
    "qsrf": {
        "units": "kg/kg",
        "long_name": "Surface Water Vapor Mass Mixing ratio (extrapolated to z=0)",
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
    "z": {"units": "m", "long_name": "full level height"},
    "pres": {"units": "Pa", "long_name": "full level pressure"},
    "presh": {"units": "Pa", "long_name": "half level pressure"},
    "u": {"units": "m/s", "long_name": "zonal wind (domain averaged)"},
    "v": {"units": "m/s", "long_name": "meridional wind (domain averaged)"},
    "t": {"units": "K", "long_name": "temperature (domain averaged)"},
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
    "omega": {
        "units": "Pa/s",
        "long_name": "large-scale pressure velocity (domain averaged)",
    },
    "T": {"units": "K", "long_name": "Liquid Water Temperature (Tl = t - (L/Cp)*ql)"},
    "q": {
        "units": "kg/kg",
        "long_name": "Total Water Mass Mixing Ratio (Vapor + Cloud Liquid)",
    },
    # need to change things here to use eulerian setup
    "divT": {
        "units": "K/s",
        "long_name": "tendency in temperature due to "
        "large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using domain mean wind ",
    },
    "xgradT": {
        "units": "K/m",
        "long_name": "gradient in temperature in zonal direction",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: used to compute eulerian horizontal advection",
    },
    "ygradT": {
        "units": "K/m",
        "long_name": "gradient in temperature in meridional direction",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: used to compute eulerian horizontal advection",
    },
    # need to change things here to use eulerian setup
    # the divq and divT terms will be computed from the x-, y-gradient instead.
    "divq": {
        "units": "kg/kg/s",
        "long_name": "tendency in water vapor due to large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using domain mean wind",
    },
    "xgradq": {
        "units": "kg/kg/m",
        "long_name": "gradient in water vapor in zonal direction",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: used to compute eulerian horizontal advection",
    },
    "ygradq": {
        "units": "kg/kg/m",
        "long_name": "gradient in water vapor in meridional direction",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: used to compute eulerian horizontal advection",
    },
    "Tref": {
        "units": "K",
        "long_name": "Reference Absolute Temperature (Mean over trajectory)",
    },
    "qref": {
        "units": "kg/kg",
        "long_name": "Reference Water Vapor Mass Mixing Ratio (Mean over trajectory)",
    },
    "qv": {"units": "kg/kg", "long_name": "water vapor mixing ratio (domain averaged)"},
    "ql": {
        "units": "kg/kg",
        "long_name": "liquid water mixing ratio (domain averaged)",
    },
    "o3mmr": {
        "units": "kg/kg",
        "long_name": "ozone mass mixing ratio (domain averaged)",
    },
    "fsdtERA": {
        "units": "W/m2",
        "long_name": "Mean top downward short-wave radiation flux",
    },
    "fsntERA": {
        "units": "W/m2",
        "long_name": "Mean top net short-wave radiation flux",
    },
    "fsntcERA": {
        "units": "W/m2",
        "long_name": "Mean top net short-wave radiation flux, clear sky",
    },
    "flntERA": {"units": "W/m2", "long_name": "Mean top net long-wave radiation flux"},
    "flntcERA": {
        "units": "W/m2",
        "long_name": "Mean top net long-wave radiation flux, clear sky",
    },
    "fsnsERA": {
        "units": "W/m2",
        "long_name": "Mean surface net short-wave radiation flux",
    },
    "fsnscERA": {
        "units": "W/m2",
        "long_name": "Mean surface net short-wave radiation flux, clear sky",
    },
    "flnsERA": {
        "units": "W/m2",
        "long_name": "Mean surface net long-wave radiation flux",
    },
    "flnscERA": {
        "units": "W/m2",
        "long_name": "Mean surface net long-wave radiation flux, clear sky",
    },
    "cloudERA": {"units": "0-1", "long_name": "cloud fraction (domain averaged)"},
    "tsec": {"units": "s", "long_name": "Time in seconds after 00Z on nbdate"},
    "calday": {"units": "d", "long_name": "'Time in days after 00Z on 31 Dec iyear"},
    "year": {"units": "year", "long_name": "Year"},
    "month": {"units": "month", "long_name": "Month"},
    "day": {"units": "day", "long_name": "Day"},
    "hour": {"units": "hour", "long_name": "Hour"},
    "nbdate": {
        "units": "yymmdd",
        "long_name": "Base date (Note that only two digit year is permitted)",
    },
    "bdate": {
        "units": "yymmdd",
        "long_name": "Base date (Note that only two digit year is permitted)",
    },
    "phis": {"units": "m2/s2", "long_name": "Surface geopotential"},
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
}


# sam variable : era5 variable
# the number of variables here need to match those in the sam_attributes;
sam_from_era5_variables = {
    "z": "height_h_local",
    "Ps": "sp_mean",
    "pres": "p_h_mean",
    "presh": "p_h_mean",
    "u": "u_mean",
    "v": "v_mean",
    "t": "t_mean",
    "qv": "q_mean",
    "ql": "clwc_mean",
    "cloudERA": "cc_mean",
    "omega": "w_pressure_corr_mean",
    "o3mmr": "o3_mean",
    "divT": "dtdt_adv",
    "divq": "dqdt_adv",
    "xgradT": "dtdx",
    "ygradT": "dtdy",
    "xgradq": "dqdx",
    "ygradq": "dqdy",
    "ug": "u_g",
    "vg": "v_g",
    "lat_traj": "lat",
    "lon_traj": "lon",
    "u_traj": "u_traj",
    "v_traj": "v_traj",
}

# era5 units : sam units
# we replace era5 units here
# added two units for the spatial gradients
era5_to_sam_units = {
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
    "K m**-1": "K/m",
    "kg kg**-1 m**-1": "kg/kg/m",
}


def from_era5(ds_era5, da_levels, parameters, metadata):
    """Obtain a sam input file from era5 variable set at high resolution"""
    # Put full levels midway between half-levels, I think this is consistent with DALES
    # Reverse order of data, to confirm to other kpt input
    # Does it apply to sam as well?
    sam_half_level_array = da_levels.values
    sam_full_level_array = 0.5 * (sam_half_level_array[:-1] + sam_half_level_array[1:])
    sam_full_level_coord = {
        "nlev": (
            "nlev",
            (np.arange(len(sam_full_level_array)) + 1.0)[::-1],
            {"long_name": "model full levels"},
        )
    }

    # note: the sam_full_level is actully height in meters, so later on in the code, I should actually
    # update the values to make it truely the pressure coord.
    sam_mean_preslev_coord = {
        "lev": (
            "lev",
            sam_full_level_array.astype(np.double),
            {"units": "Pa", "long_name": "Pressure level (time-averaged)"},
        )
    }
    # nDS_coord = {"nDS": ("nDS", [0], {},)}
    # sam needs lat and lon as well in the coordinate
    sam_lon_coord = {
        "lon": (
            "lon",
            [ds_era5.lon.values[0].astype(np.double)],
            {},
        )
    }
    sam_lat_coord = {
        "lat": (
            "lat",
            [ds_era5.lat.values[0].astype(np.double)],
            {},
        )
    }

    ds_sam = xr.Dataset(
        coords={
            "time": ds_era5.time,
            **sam_full_level_coord,
            **sam_mean_preslev_coord,
            **sam_lat_coord,
            **sam_lon_coord,
        }
    )
    # Variables from dictionary
    # Including unit checks
    for variable in sam_from_era5_variables:
        era5_var = sam_from_era5_variables[variable]
        da_era5 = ds_era5[era5_var]
        # perform units check
        unit_guess = era5_to_sam_units.get(da_era5.units, da_era5.units)
        if not unit_guess == sam_attributes[variable]["units"]:
            except_str = (
                "Incompatible units between ERA5 and sam for variable "
                + variable
                + ". Please fix using the fix_era5_to_sam_units "
                + "dictionary: ERA converted variable is "
                + unit_guess
                + ", sam variable is "
                + sam_attributes[variable]["units"]
            )
            raise Exception(except_str)
        # single level variable
        if np.ndim(da_era5.values) == 1:
            (Ntime,) = da_era5.values.shape
            # prescribe the dimensions of the variable as well as following:
            ds_sam[variable] = (
                ("time", "lat", "lon"),
                da_era5.values.reshape((Ntime, 1, 1)).astype(np.single),
                sam_attributes[variable],
            )
        # half level variable
        elif variable in ["zh", "presh"]:
            da_era5_on_half_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, sam_half_level_array
            )
            Ntime, Nlev1 = da_era5_on_half_levels.shape
            ds_sam[variable] = (
                ("time", "nlevp1", "lat", "lon"),
                da_era5_on_half_levels.reshape((Ntime, Nlev1, 1, 1)).astype(np.single),
                sam_attributes[variable],
            )
        # full level variable
        else:
            da_era5_on_full_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, sam_full_level_array
            )
            #  get dimensions out:
            Ntime, Nlev = da_era5_on_full_levels.shape
            ds_sam[variable] = (
                ("time", "lev", "lat", "lon"),
                da_era5_on_full_levels.reshape((Ntime, Nlev, 1, 1)).astype(np.single),
                sam_attributes[variable],
            )

    variables_to_centralise = {
        "fsnsERA": "msnswrf_mean",
        "flnsERA": "msnlwrf_mean",
        "fsntERA": "mtnswrf_mean",
        "flntERA": "mtnlwrf_mean",
        "fsntcERA": "mtnswrfcs_mean",
        "flntcERA": "mtnlwrfcs_mean",
        "fsnscERA": "msnswrfcs_mean",
        "flnscERA": "msnlwrfcs_mean",
        "fsdtERA": "mtdwswrf_mean",
    }
    for variable in variables_to_centralise:
        this_central_estimate = central_estimate(
            ds_era5[variables_to_centralise[variable]].values
        )
        this_central_estimate = this_central_estimate.reshape(
            (
                Ntime,
                1,
                1,
            )
        )
        ds_sam[variable] = (
            ("time", "lat", "lon"),
            this_central_estimate.astype(np.single),
            sam_attributes[variable],
        )

    # Manually compute the needed variables for SAMIOP, which are not directly available from ERA5 #
    # Tl, qtot:
    # J/kg specific heat of condensation
    L = hrv - hrl

    T_lw = ds_sam["t"].values - L / cpd * ds_sam["ql"].values
    ds_sam["T"] = (
        ("time", "lev", "lat", "lon"),
        T_lw.astype(np.single),
        sam_attributes["T"],
    )

    qtot = ds_sam["qv"].values + ds_sam["ql"].values
    ds_sam["q"] = (
        ("time", "lev", "lat", "lon"),
        qtot.astype(np.single),
        sam_attributes["q"],
    )

    # Eulerian version of horizontal T, q advection:
    divT = -(ds_sam["u"] * ds_sam["xgradT"] + ds_sam["v"] * ds_sam["ygradT"]).values
    divq = -(ds_sam["u"] * ds_sam["xgradq"] + ds_sam["v"] * ds_sam["ygradq"]).values

    ds_sam["divT"] = (
        ("time", "lev", "lat", "lon"),
        divT.astype(np.single),
        sam_attributes["divT"],
    )

    ds_sam["divq"] = (
        ("time", "lev", "lat", "lon"),
        divq.astype(np.single),
        sam_attributes["divq"],
    )

    # surface pressure: (from ERA5 surface pressure)
    Ps = ds_sam["presh"][:, -1, :, :].values
    # ds_sam["Ps"] = (
    #    ("time", "lat", "lon"),
    #    Ps.astype(np.single).reshape((Ntime, 1, 1)),
    #    sam_attributes["Ps"],
    # )
    # surface temperature and humidity: (needs interpolation)
    Nt = len(ds_era5.time)
    T_surf = np.zeros((Nt, 1, 1))
    q_surf = np.zeros((Nt, 1, 1))
    for it in range(Nt):
        FT = intp.interp1d(
            ds_sam["pres"][it, :, 0, 0].values,
            ds_sam["t"][it, :, 0, 0].values,
            bounds_error=False,
            fill_value="extrapolate",
        )
        Fq = intp.interp1d(
            ds_sam["pres"][it, :, 0, 0].values,
            ds_sam["qv"][it, :, 0, 0].values,
            bounds_error=False,
            fill_value="extrapolate",
        )
        T_surf[it, :, :] = FT(Ps[it, 0, 0])
        q_surf[it, :, :] = Fq(Ps[it, 0, 0])

    ds_sam["Tsair"] = (
        ("time", "lat", "lon"),
        T_surf.astype(np.single),
        sam_attributes["Tsair"],
    )
    ds_sam["qsrf"] = (
        ("time", "lat", "lon"),
        q_surf.astype(np.single),
        sam_attributes["qsrf"],
    )

    # remove surface omega (based on the notion that surface vertical velocity corresponds to the atmos. tide)
    omega = ds_sam["omega"][:, :, 0, 0].values
    omega_out = np.zeros_like(omega)
    # time averaged pressure levels
    lev = ds_sam["pres"][:, :, 0, 0].values.mean(axis=0)
    for it in range(Nt):
        xx = np.maximum(
            np.zeros_like(lev),
            np.minimum(np.ones_like(lev), (700e2 - lev) / (700e2 - 150e2)),
        )
        f = 0.5 * (1 + np.cos(np.pi * xx))
        omega_out[it, :] = omega[it, :] - f * omega[it, -1]

    Nt, Nlev = omega_out.shape
    ds_sam["omega"] = (
        ("time", "lev", "lat", "lon"),
        omega_out.astype(np.single).reshape((Nt, Nlev, 1, 1)),
        sam_attributes["omega"],
    )

    # update the pressure coordinate values
    ds_sam["lev"] = (
        ("lev"),
        lev.astype(np.double),
        {"units": "Pa", "long_name": "Pressure level (time-averaged)"},
    )

    # surface pressure tendencey (zeros, equiv. to surface omega, which has been removed above.)
    ds_sam["Ptend"] = (
        ("time", "lat", "lon"),
        np.zeros_like(ds_sam["qsrf"].values).astype(np.single),
        sam_attributes["Ptend"],
    )

    # build Tref and qref: time or trajectory averaged t and q
    Nt, nlev = np.shape(ds_sam["t"][:, :, 0, 0])

    tvec = np.ones((Nt, 1))
    Tref_vec = np.zeros((1, nlev))
    qref_vec = np.zeros((1, nlev))

    Tref_vec[0, :] = ds_sam["t"][:, :, 0, 0].mean(axis=0).values
    Tref = Tref_vec * tvec

    qref_vec = ds_sam["qv"][:, :, 0, 0].mean(axis=0).values
    qref = qref_vec * tvec

    ds_sam["Tref"] = (
        ("time", "lev", "lat", "lon"),
        Tref.reshape((Nt, nlev, 1, 1)).astype(np.single),
        sam_attributes["Tref"],
    )
    ds_sam["qref"] = (
        ("time", "lev", "lat", "lon"),
        qref.astype(np.single).reshape((Nt, nlev, 1, 1)),
        sam_attributes["qref"],
    )
    Tg_skin = ds_era5["skt_mean"].values
    ds_sam["Tg"] = (
        ("time", "lat", "lon"),
        Tg_skin.reshape((Nt, 1, 1)).astype(np.single),
        sam_attributes["Tg"],
    )
    # Surface fluxes: obtain from time mean in ERA data, change sign for SAM
    sfc_sens_flx = central_estimate(ds_era5["msshf_mean"].values)
    ds_sam["shflx"] = (
        ("time", "lat", "lon"),
        -sfc_sens_flx.reshape((Nt, 1, 1)).astype(np.single),
        sam_attributes["shflx"],
    )
    sfc_lat_flx = central_estimate(ds_era5["mslhf_mean"].values)
    ds_sam["lhflx"] = (
        ("time", "lat", "lon"),
        -sfc_lat_flx.reshape((Nt, 1, 1)).astype(np.single),
        sam_attributes["lhflx"],
    )
    # Final checks: are all variables present?
    ds_sam["time_traj"] = (
        ds_era5["time"] - np.datetime64("1970-01-01T00:00")
    ) / np.timedelta64(1, "s")
    ds_sam["time_traj"].attrs.update(**sam_attributes["time_traj"])
    ds_sam["DS"] = (("nDS"), ["Trajectory origin"], sam_attributes["DS"])
    ds_sam["timDS"] = (
        ("nDS"),
        [
            (ds_era5["origin_datetime"] - np.datetime64("1970-01-01T00:00"))
            / np.timedelta64(1, "s")
        ],
        sam_attributes["timDS"],
    )
    ds_sam["latDS"] = (("nDS"), [ds_era5["origin_lat"]], sam_attributes["latDS"])
    ds_sam["lonDS"] = (("nDS"), [ds_era5["origin_lon"]], sam_attributes["lonDS"])

    # Follow Peter Blossey's netCDF formatting.
    # add extra dimenson variables: (time, calday, year, month, day, hour, nbdate, bdate, phis)
    # note: these variables already satisfied the data type

    # set up the right value for "tsec" dimension according to Peter Blossey's script
    # --------------------------------------------------------------------------- #
    sam_timevars = convert_era5_time_to_sam_timevars(ds_sam["time"])
    # --------------------------------------------------------------------------- #
    Dimension_variables = {
        "tsec": [np.int32, ("time"), sam_timevars["tsec"]],
        "calday": [np.double, ("time"), sam_timevars["calday"]],
        "year": [np.int32, ("time"), sam_timevars["year"]],
        "month": [np.int32, ("time"), sam_timevars["month"]],
        "day": [np.int32, ("time"), sam_timevars["day"]],
        "hour": [np.double, ("time"), sam_timevars["hour"]],
        "nbdate": [np.int32, (), sam_timevars["nbdate"]],
        "bdate": [np.int32, (), sam_timevars["nbdate"]],
        "phis": [np.double, ("lat", "lon"), np.zeros((1, 1)).astype(np.double)],
    }

    for variable in Dimension_variables:
        var_dtype = Dimension_variables[variable][0]
        var_dim = Dimension_variables[variable][1]
        var_val = Dimension_variables[variable][2]
        ds_sam[variable] = (
            var_dim,
            var_val.astype(var_dtype),
            sam_attributes[variable],
        )
        #   update the attribute for calday:
        if variable == "calday":
            yyyy = Dimension_variables["year"][2]
            ds_sam["calday"] = ds_sam["calday"].assign_attrs(
                {
                    "long_name": "Time in days after 00Z on 31 Dec {0:s}".format(
                        str(yyyy[0] - 1)
                    )
                }
            )

    # Change order of data, to confirm to other sam input
    ds_sam = ds_sam.sortby("lev", ascending=True)
    ds_sam = ds_sam.sortby("nlevp1", ascending=True)
    for var in sam_attributes:
        if var not in ds_sam:
            print(var + " is missing in the sam formatted output")
    # Needs improvement still
    sam_dict = {
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
        # "t_skin_correct": "Skin temperature has been corrected "
        # "by 1.000000. Motivation: value from IFS is actually "
        # "the open SST, which is lower than the skin temperature.",
        "omega_correct": "surface omega corresponds to the atmospheric tide has been removed."
        "by assuming its vertical structure is uniformly one up to 700 hPa and zero above 150hPa (P.Blossey)",
    }
    ds_sam.attrs.update(**sam_dict)
    return ds_sam


def convert_era5_time_to_sam_timevars(era5_time):
    """
    Note: era5_time is a xarray data array
          calday: calendar day (or day or years) in SAM (1~365 or 366), double
          tsec: seconds from the base date "nbdate" at 0Z. this is second of day if
                simulation only last within 1 day.
    """
    refyear = (era5_time.dt.strftime("%Y")).values.astype(np.int32)[0]
    refdate = np.datetime64(str(refyear - 1) + "-12-31")
    calday = (era5_time.values - refdate) / np.timedelta64(86400, "s")

    # time in second from the base date "nbdate"
    tsec = np.round(86400 * (calday - np.floor(calday[0])))

    return {
        "calday": calday.astype(np.double),
        "tsec": tsec.astype(np.int32),
        "year": (era5_time.dt.strftime("%Y")).values.astype(np.int32),
        "month": (era5_time.dt.strftime("%m")).values.astype(np.int32),
        "day": (era5_time.dt.strftime("%d")).values.astype(np.int32),
        "hour": (era5_time.dt.strftime("%H")).values.astype(np.double),
        "nbdate": (era5_time.dt.strftime("%y%m%d")).values.astype(np.int32)[0],
    }
