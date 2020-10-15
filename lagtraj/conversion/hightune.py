"""Module describing the hightune variables as a dictionary"""

hightune_variables = {
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
