"""Module describing the RACMO variables as a dictionary"""

racmo_variables = {
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
        "t_skin_correct": "Skin temperature has been corrected by 1.000000. "
        "Motivation: value from IFS is actually the open SST, "
        "which is lower than the skin temperature.",
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
        "units": "0-1",
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
}
