"""This module contains:
    1. a list of the SAM variables as a dictionary
    2. the mapping between the internal lagtraj variables (era5) and SAM
    3. and a function to process the conversion
To-DO: correct the variable dimension to match the assigned one. 
       broadcast to have [time, nlev, 1, 1] 
"""
import xarray as xr
import numpy as np
import scipy.interpolate as intp
import datetime


from ....domain.sources.era5.constants import rg
from ....utils.interpolation.methods import (
    steffen_1d_no_ep_time,
    central_estimate,
)

# refer to Peter's matlab script; It seems that it works well with the kpt formatted output as well.
sam_attributes = {
## following Peter Blossey's matlab script structure:
# 1. Surface level variables (manually added)
    "Ps": {"units": "Pa", "long_name": "surface pressure"},
    "Ptend" : {"units": "Pa/s", "long_name":"surface pressure tendency"},
    "Tg": {
        "units": "K",
        "long_name": "Surface (skin) Temperature (SST if over water)",
        "t_skin_correct": "Skin temperature has been corrected by 1.000000. "
        "Motivation: value from IFS is actually the open SST, "
        "which is lower than the skin temperature.",
    },
    # both value needs to be change sign (done)
    "shflx": {"units": "W/m2", "long_name": "surface sensible heat flux"},
    "lhflx": {"units": "W/m2", "long_name": "surface latent heat flux"},
    # thest two needs to be computed  (done)
    "Tsair" :{"units": "K", "long_name": "Surface Air Temperature (extrapolated to z=0)"},
    "qsrf" : {"units": "kg/kg", "long_name": "Surface Water Vapor Mass Mixing ratio (extrapolated to z=0)"},
    # need to double check which is the right lon, lat for SAM and how they are different.
    "lat": {"units": "degrees North", "long_name": "latitude"},
    "lon": {"units": "degrees East", "long_name": "longitude"},
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
# 2. Full level (vertical varying) variables:
    "zf": {"units": "m", "long_name": "full level height"},
#    "zh": {"units": "m", "long_name": "half level height"},
    "pres": {"units": "Pa", "long_name": "full level pressure"},
    "presh": {"units": "Pa", "long_name": "half level pressure"},
    "u": {"units": "m/s", "long_name": "zonal wind (domain averaged)"},
    "v": {"units": "m/s", "long_name": "meridional wind (domain averaged)"},
#    "t": {"units": "K", "long_name": "temperature (domain averaged)"},
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
# Note: these two variables need to be calculated from the era5 variables.
    "T" : {"units": "K", "long_name": "Liquid Water Temperature (Tl = T - (L/Cp)*ql)"}, #"standard_name": "air_temperature"}, # Tl
    "q" : {"units":"kg/kg", "long_name": "Total Water Mass Mixing Ratio (Vapor + Cloud Liquid)"}, # qtot
    "divT": {
        "units": "K/s",
        "long_name": "tendency in temperature due to "
        "large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using velocity relative to wind "
        "on trajectory (u_traj,v_traj)",
    },
    "divq": {
        "units": "kg/kg/s",
        "long_name": "tendency in water vapor due to large-scale horizontal advection",
        "info": "derived at pressure levels",
        "lagrangian": "Lagrangian setup: horizontal advection calculated "
        "using velocity relative to wind on trajectory (u_traj,v_traj)",
    },
    "Tref": {"units": "K", "long_name": "Reference Absolute Temperature (Mean over trajectory)"},
    "qref": {"units": "kg/kg", "long_name": "Reference Water Vapor Mass Mixing Ratio (Mean over trajectory)"},
    "qv": {"units": "kg/kg", "long_name": "water vapor mixing ratio (domain averaged)"},
    "ql": {
        "units": "kg/kg",
        "long_name": "liquid water mixing ratio (domain averaged)",
    },
#    "qi": {"units": "kg/kg", "long_name": "ice water mixing ratio (domain averaged)"},
    "o3mmr": {"units": "kg/kg", "long_name": "ozone mass mixing ratio (domain averaged)"},
  # 3. ERA radiative flux outputs:
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
# 4. cloud fraction
    "cloudERA": {"units": "0-1", "long_name": "cloud fraction (domain averaged)"},
# 5.dimension related parameters:
    "lev": {"units": "Pa", "long_name": "Pressure level"},
    "tsec": {"units": "s", "long_name": "Time in seconds after 00Z on nbdate"},
    "time": {"units": "s", "long_name": "Time in seconds after 00Z on nbdate"},
#    "calday": {"units": "d", "long_name": "'Time in days after 00Z on 31 Dec iyear"},
    "year": {"units": "year", "long_name": "Year"},
    "month": {"units": "month", "long_name": "Month"},
    "day": {"units": "day", "long_name": "Day"},
    "hour": {"units": "hour", "long_name": "Hour"},
    "nbdate": {"units": "yymmdd", "long_name": "Base date (Note that only two digit year is permitted)"},
    "bdate": {"units": "yymmdd", "long_name": "Base date (Note that only two digit year is permitted)"},
    "phis": {"units": "m2/s2", "long_name": "Surface geopotential"},
}


# sam variable : era5 variable
# (we loop over sam variables here)
# XYC Development Note: I think the first column is the variable name in SAM and the second is the 
# era5 variable name
# the number of variables here need to match those in the sam_attributes;
sam_from_era5_variables = {
    "zf": "height_h_local",
#    "zh": "height_h_local",
#    "ps": "sp_mean",
    "pres": "p_h_mean",
    "presh": "p_h_mean",
    "u": "u_mean",
    "v": "v_mean",
#    "t": "t_mean",
    "qv": "q_mean",
    "ql": "clwc_mean",
#    "qi": "ciwc_mean",
    "cloudERA": "cc_mean",
    "omega": "w_pressure_corr_mean",
    "o3mmr": "o3_mean",
#    "t_local": "t_local",
#    "q_local": "q_local",
#    "ql_local": "clwc_local",
#    "qi_local": "ciwc_local",
#    "u_local": "u_local",
#    "v_local": "v_local",
#    "cc_local": "cc_local",
    "tadv": "dtdt_adv",
    "qadv": "dqdt_adv",
#    "uadv": "dudt_adv",
#    "vadv": "dvdt_adv",
    "ug": "u_g",
    "vg": "v_g",
#    "tladv": "dt_ldt_adv",
#    "qladv": "dclwcdt_adv",
#    "qiadv": "dciwcdt_adv",
#    "ccadv": "dccdt_adv",
    "lat": "lat",
    "lon": "lon",
    "lat_traj": "lat",
    "lon_traj": "lon",
#    "lat_grid": "lat",
#    "lon_grid": "lon",
    "u_traj": "u_traj",
    "v_traj": "v_traj",
#    "open_sst": "sst_mean",
}

# era5 units : sam units
# we replace era5 units here
era5_to_sam_units = {
    "m s**-1": "m/s",
    "1": "0-1",
    "degrees_north": "degrees North",
    "degrees_east": "degrees East",
    "metres": "m",
    "m": "m",
    "kg kg**-1": "kg/kg",
    "Pa s**-1": "Pa/s",
    "K s**-1": "K/s",
    "kg kg**-1 s**-1": "kg/kg/s",
    "m s**-1 s**-1": "m/s2",
    "kg m**-3": "kg/m3",
    "s**-1": "1/s",
    "W m**-2": "W/m2",
}


def from_era5(ds_era5, da_levels, parameters, metadata):
    """Obtain a sam input file from era5 variable set at high resolution"""
    # Put full levels midway between half-levels, I think this is consistent with DALES
    # Reverse order of data, to confirm to other kpt input
    ### Does it apply to sam as well? 
    sam_half_level_array = da_levels.values
    sam_full_level_array = 0.5 * (sam_half_level_array[:-1] + sam_half_level_array[1:])
#    sam_half_level_coord = {
#        "nlevp1": (
#            "nlevp1",
#            (np.arange(len(sam_half_level_array)) + 1.0)[::-1],
#            {"long_name": "model half levels"},
#        )
#    }
    sam_full_level_coord = {
        "nlev": (
            "nlev",
            (np.arange(len(sam_full_level_array)) + 1.0)[::-1],
            {"long_name": "model full levels"},
        )
    }
    sam_mean_preslev_coord = {
        "lev": (
            "lev", 
            ds_era5["p_h_mean"].values.mean(axis=0),  # time averaged pressure levels
            {"units":"Pa", "long_name": "Pressure level (time-averaged)"} ) 
    }
    nDS_coord = {"nDS": ("nDS", [0], {},)}
#    sam_soil_coord = {
#        "nlevs": ("nlevs", np.arange(4) + 1.0, {"long_name": "soil levels"})
#    }
#   sam needs lat and lon as well in the coordinate
    sam_lon_coord = {
        "lon": (
            "lon", 
            ds_era5.lon[0].values,
            {"units": "degrees_north", "long_name": "Longitude"},
        )
    }
    sam_lat_coord = {
        "lat": (
            "lat", 
            ds_era5.lat[0].values,
            {"units": "degrees_east", "long_name": "Latitude"},
        )
    }

    # set up the right value for "tsec" dimension according to Peter Blossey's script
    # this part here can be turned in to a function in this python script. 
    # --------------------------------------------------------------------------- #
    sam_timevars = convert_era5_time_to_sam_timevars(ds_era5.time.values)
    sam_time_coord = {
        "tsec": (
            "tsec", 
            sam_timevars['time'],
            {"units": "s", "long_name": "Time in seconds after 00Z on nbdate"},
        )
    } 
    # --------------------------------------------------------------------------- #

    ds_sam = xr.Dataset(
        coords={
#            "tsec": sam_timevars['time'], #ds_era5.time,
            **sam_time_coord,
            **sam_full_level_coord,
            **sam_mean_preslev_coord,
#            **sam_half_level_coord,
            **sam_lat_coord,
            **sam_lon_coord,
#            **sam_soil_coord,
#            **nDS_coord,
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
            Ntime = da_era5.values.shape
            # prescribe the dimensions of the variable as well as following:
            ds_sam[variable] = (
                ("tsec","lat", "lon"), 
                da_era5.values.reshape((Ntime, 1,1)).astype(np.single),
                sam_attributes[variable],
            ) 
        # half level variable
        elif variable in ["zh", "presh"]:
            da_era5_on_half_levels = steffen_1d_no_ep_time(
                da_era5.values, ds_era5["level"].values, sam_half_level_array
            )
            Ntime, Nlev1 = da_era5_on_half_levels.shape
            ds_sam[variable] = (
                ("tsec", "nlevp1","lat","lon"),
                da_era5_on_half_levels.reshape((Ntime, Nlev1, 1,1)).astype(np.single),
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
                ("tsec", "lev", "lat", "lon"),
                da_era5_on_full_levels.reshape((Ntime, Nlev, 1, 1)).astype(np.single),
                sam_attributes[variable],
            )

    # Simple unit fix fails for these variables
    # So these are added manually after checking
    # that units are compatible
#    variables_to_manually_add = {
#        "high_veg_type": "tvh_mean",
#        "low_veg_type": "tvl_mean",
#        "high_veg_lai": "lai_hv_mean",
#        "low_veg_lai": "lai_lv_mean",
#        "slor": "slor_mean",
#        "q_skin": "src_mean",
#        "snow": "sd_mean",
#        "lsm": "lsm_mean",
#    }
#    for variable in variables_to_manually_add:
#        ds_sam[variable] = ds_era5[variables_to_manually_add[variable]]
#        ds_sam[variable].attrs.update(**sam_attributes[variable])
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
        this_central_estimate = this_central_estimate.reshape((Ntime,1,1,))
        ds_sam[variable] = (
            ("tsec","lat","lon"),
            this_central_estimate.astype(np.single),
            sam_attributes[variable],
        )
#### XYC note: I think this is where I can compute extra varialbles from ERA5 and add to ds_sam
    # Tl, qtot:
    L = 2.5e6          # J/kg specific heat of condensation
    Cp = 1004          # specific heat of air at constant pressure p; J/kg/Kcontanst as in SAM
    T_lw = ds_era5['t_mean'].values - L/Cp * ds_sam['ql'].values
    Ntime, Nlev = T_lw.shape
    ds_sam["T"] = (
        ("tsec", "lev","lat","lon"), 
        T_lw.astype(np.single).reshape((Ntime,Nlev,1,1)),
        sam_attributes["T"],
    )
        
    qtot = ds_sam['qv'].values + ds_sam['ql'].values
    ds_sam["q"] = (
        ("tsec", "lev","lat","lon"), 
        qtot.astype(np.single).reshape((Ntime, Nlev, 1,1)),
        sam_attributes["q"],
    )

    # surface pressure:
    Ps = ds_sam["presh"][:,-1].values
    ds_sam['Ps'] = (
        ("tsec","lat","lon"), 
        Ps.astype(np.single).reshape((Ntime, 1,1)),  
        sam_attributes["Ps"], 
    )
    # surface temperature and humidity: (needs interpolation)
    Nt = len(ds_era5.time)
    T_surf = np.zeros((Nt,1,1))
    q_surf = np.zeros((Nt,1,1))
    for it in range(Nt):    
        FT = intp.interp1d(ds_sam['pres'][it,:].values, ds_sam['t'][it,:].values, bounds_error=False, fill_value='extrapolate')
        Fq = intp.interp1d(ds_sam['pres'][it,:].values, ds_sam['qv'][it,:].values, bounds_error=False, fill_value='extrapolate')
        T_surf[it,:,:] = FT(Ps[it])
        q_surf[it,:,:] = Fq(Ps[it])
    
    ds_sam["Tsair"] = (
        ("tsec","lat","lon"), 
        T_surf.astype(np.single),
        sam_attributes["Tsair"],
    )
    ds_sam["qsrf"] = (
        ("tsec","lat","lon"), 
        q_surf.astype(np.single),
        sam_attributes["qsrf"],
    )

    # remove surface omega (based on the notion that surface vertical velocity corresponds to the atmos. tide)
    omega = ds_sam["omega"].values
    omega_out = np.zeros(np.shape(omega))
    lev = ds_sam["pres"].values.mean(axis=0)  # time averaged pressure levels
    for it in range(Nt):
        xx = np.maximum(np.zeros_like(lev), np.minimum(np.ones_like(lev), (700e2-lev)/(700e2-150e2) ))
        f = 0.5*(1 + np.cos(np.pi * xx))  # shape (nlev,)
        omega_out[it,:] = omega[it,:] - f*omega[it,-1]

    Ntime, Nlev = omega_out.shape
    ds_sam["omega"] = (
        ("tsec", "lev","lat","lon"), 
        omega_out.astype(np.single).reshape((Ntime, Nlev, 1,1)),
        sam_attributes["omega"],
    )

    # surface pressure tendencey (zeros, equiv. to surface omega, which has been removed above.)
    ds_sam["Ptend"] = (
        ("tsec", "lev","lat","lon"), 
        np.zeros_like(omega_out).astype(np.single), 
        sam_attributes["Ptend"]
    )
    
    # build Tref and qref: time or trajectory averaged t and q
    Nt, nlev = np.shape(ds_era5['t_mean'])
    tvec = np.ones((Nt, 1))
    Tref_vec = np.zeros((1, nlev))
    qref_vec = np.zeros((1, nlev))

    Tref_vec[0,:] = ds_era5['t_mean'].mean(axis=0).values #(broadcast to have shape time, nlev)
    Tref = Tref_vec * tvec
    
    qref_vec[0,:] = ds_sam['qv'].mean(axis=0).values
    qref = qref_vec * tvec

    ds_sam["Tref"] = (
        ("tsec", "lev", "lat", "lon"), 
        Tref.reshape((Nt, nlev, 1,1)).astype(np.single), 
        sam_attributes["Tref"]
    )
    ds_sam["qref"] = (
        ("tsec", "lev", "lat", "lon"), 
        qref.reshape((Nt, nlev, 1,1)).astype(np.single), 
        sam_attributes["qref"]
    )

    # Soil moisture: combine levels
#    swvl1 = ds_era5["swvl1_mean"].values
#    swvl2 = ds_era5["swvl2_mean"].values
#    swvl3 = ds_era5["swvl3_mean"].values
#    swvl4 = ds_era5["swvl4_mean"].values
#    q_soil = np.stack((swvl1, swvl2, swvl3, swvl4), axis=-1)
#    ds_sam["q_soil"] = (
#        ("time", "nlevs"),
#        q_soil,
#        sam_attributes["q_soil"],
#    )
    # Soil temperature: combine levels
#    stl1 = ds_era5["stl1_mean"].values
#    stl2 = ds_era5["stl2_mean"].values
#    stl3 = ds_era5["stl3_mean"].values
#    stl4 = ds_era5["stl4_mean"].values
#    t_soil = np.stack((stl1, stl2, stl3, stl4), axis=-1)
#    ds_sam["t_soil"] = (
#        ("time", "nlevs"),
#        t_soil,
#        sam_attributes["t_soil"],
#    )
    # Soil thickness: combine levels
#    h_soil = np.array([0.07, 0.21, 0.72, 1.89])
#    ds_sam["h_soil"] = (
#        ("nlevs"),
#        h_soil,
#        sam_attributes["h_soil"],
#    )
    # Orography: derive from surface geopotential
#    ds_sam["orog"] = ds_era5["z_mean"] / rg
#    ds_sam["orog"].attrs.update(**sam_attributes["orog"])
    # Heat roughness, derive from "flsr" variable
#    ds_sam["heat_rough"] = np.exp(ds_era5["flsr_mean"])
#    ds_sam["heat_rough"].attrs.update(**sam_attributes["heat_rough"])

    # Apply correction to t_skin (see output files)
#    ds_sam["Tg"] = ds_era5["stl1_mean"].values + 1.0
#    ds_sam["Tg"].attrs.update(**sam_attributes["Tg"])
    Tg_corr = ds_era5["stl1_mean"].values + 1.0
    ds_sam["Tg"] = (
        ("tsec", "lat","lon"), 
        Tg_corr.reshape((Nt, 1,1)).astype(np.single),
        sam_attributes["Tg"],
    )
    # Surface fluxes: obtain from time mean in ERA data, do not change sign (do change sign for SAM)
    sfc_sens_flx = central_estimate(ds_era5["msshf_mean"].values)
    ds_sam["shflx"] = (
        ("tsec", "lat", "lon"),
        -sfc_sens_flx.reshape((Nt, 1,1)).astype(np.single),
        sam_attributes["shflx"],
    )
    sfc_lat_flx = central_estimate(ds_era5["mslhf_mean"].values)
    ds_sam["lhflx"] = (
        ("tsec", "lat", "lon"),
        -sfc_lat_flx.reshape((Nt,1,1)).astype(np.single),
        sam_attributes["lhflx"],
    )
    # Final checks: are all variables present?
    # also putting in several time variables and lat lon variables:
    # the time structuring may need to be changed for SAM:
    # fix this tomorrow;
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
  
    # Peter Blossey' computed iyear and calday for the time variable and its attributes etc. Need to add those to here as well.
    # add dimenson variables: (lat, lon, lev, tsec, time, calday, year, month, day, hour, nbdate, bdate, phis)
     # varn: [datatype, dimensions]
    Dimension_variables = {
#        "lev": [np.double, ("nlev"), lev],
#        "tsec": [np.int32, ("tsec"), time],
        "time": [np.int32, ("tsec"), sam_timevars['time']],
        "calday": [np.double, ("tsec"), sam_timevars['calday']],
        "year": [np.int32, ("tsec"), sam_timevars['year']],
        "month": [np.int32, ("tsec"), sam_timevars['month']],
        "day": [np.int32, ("tsec"), sam_timevars['day']],
        "hour": [np.double, ("tsec"), sam_timevars['hour']],
        "nbdate": [np.int32, (), sam_timevars['nbdate']],
        "bdate": [np.int32, (), sam_timevars['nbdate']],
        "phis": [np.double, ("lat","lon"), 0],
    }
    sam_dimvars_attrs = {
        "time": {"units":"s", "long_name": "Time in seconds after 00Z on nbdate"},
        "calday": {"units": "d", "long_name": "Time in days after 00Z on 31 Dec"},
        "year" : {"units": "year", "long_name": "Year of data"},
        "month": {"units": "month", "long_name": "Month of data"},
        "day"  : {"units": "day", "long_name": "day of data"}, 
        "hour" : {"units": "hour", "long_nmae": "data hour"},
        "nbdate" : {"units": "yyyymmdd", "long_name": "Base date (Note that only two digit year is permitted)"},
        "bdate" : {"units": "yyyymmdd", "long_name": "Base date (Note that only two digit year is permitted)"},
        "phis" : {"units": "m2/s2", "long_name": "Surface geopotential"}
    }

    for variable in Dimension_variables:
        var_dtype = Dimension_variables[variable][0]
        var_dim = Dimension_variables[variable][1]
        var_val = Dimension_variables[variable][2] 
        ds_sam[variable] = (
            var_dim,
            var_val.astype(var_dtype),  # to be created in a function
            sam_dimvars_attrs[variable], # to be added
        )
        #   update the attribute for calday:
        if variable=='calday':
           yyyy = Dimension_variables['year'][2]
           ds_sam["calday"].attrs.update=({"long_name": 
                  "Time in days after 00Z on 31 Dec {0:s}".format(str(yyyy[0]-1))
           })


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
        "t_skin_correct": "Skin temperature has been corrected "
        "by 1.000000. Motivation: value from IFS is actually "
        "the open SST, which is lower than the skin temperature.",
        "omega_correct" : "surface omega corresponds to the atmospheric tide has been removed."
        "by assuming its vertical structure is uniformly one up to 700 hPa and zero above 150hPa (P.Blossey)"
    }
    ds_sam.attrs.update(**sam_dict)
    return ds_sam


def convert_era5_time_to_sam_timevars(era5_time):
    """
     Note: era5_time is a numpy data array.
    """
    yyyy = era5_time.astype('datetime64[Y]').astype(int) + 1970
    mm = era5_time.astype('datetime64[M]').astype(int) % 12 + 1
    dd = era5_time.astype('datetime64[D]') - era5_time.astype('datetime64[M]') + 1
    Hour = era5_time.astype('datetime64[h]') - era5_time.astype('datetime64[D]')
    # need to add loop, how can I do it without loop?
    calday = np.zeros(np.shape(era5_time))
    for i in range(len(yyyy)):
        refdates = np.datetime64(str(yyyy[i]-1) + '-12-31')
        calday[i] = (era5_time[i] - refdates)/np.timedelta64(86400, 's')
    time_in_sec = np.round(86400*(calday - np.floor(calday[0])))
    time_in_sec = time_in_sec.astype(np.int32)
        
    #time = np.round(86400*(calday - np.floor(calday[0]))) same as time_in_sec above.

    yy = yyyy[0]%100
    nbdate = 1e4*yy + 1e2*mm[0] + dd[0]

    return {'calday': calday, 'time': time_in_sec,'year': yyyy, 'month': mm,
            'day': dd, 'hour': Hour, 'nbdate': nbdate}
