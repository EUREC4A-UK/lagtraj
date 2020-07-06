"""
ERA5 utilities that can
- Add heights and pressures to an input data array on model levels
- Interpolate from model levels to constant height levels (using Steffen interpolation)
- Calculate gradients using boundary values or regression method
- Extract local profiles and mean profiles
- Subselect a domain
- Filter/mask data: e.g. "ocean values only"
- Add auxiliary variables

TODO
- More focings/variable (e.g. humidity forcings and large-scale tendencies).
- Use of local versus averaged velocity in trajectories
- Distinguish box sizes for averages, gradients and trajectory calculations.
- Implement selection as ciricle with radius, rather than box? M
- Optimise code (note that interpolation is currently expensive, possibly because coordinates are not assumed to be ordered)
- Add more auxiliary variables
- Move some functionality (e.g. auxiliary variables) to more generic utilities
- Test/develop way of dealing with antimeridian (and poles, for circular boxes).
  Note: HIGH-TUNE/DEPHY prefers longitude between -180..180?
- Implement HIGH-TUNE conventions, variable renaming, attributes
- HIGH-TUNE/DEPHY needs netcdf3?
- Discuss need to check/convert float vs double (HIGH-TUNE/DEPHY expects double)
- Discuss data filter and weight procedures
- Use more exact location for means and gradients based on interpolation/weights?
- Use budgets over boundarie instead of gradients for large-scale tendencies?
- Compare regression and boundary gradients for data
- Look into other mean/gradient techniques (e.g. Gaussian weighted, see CSET code).
- Further documentation
- Keep checking against cf conventions
"""

import os
import numbers
import datetime
import numpy as np
import pandas as pd
import xarray as xr


def add_auxiliary_variable(ds_to_expand, var, settings_dictionary):
    """Adds auxiliary variables to arrays.
    Alternatively, the equations could be separated out to another utility
    I think this may be adding a 'black box layer' though
    To be discussed"""
    if var == "theta":
        attr_dict = {"units": "K", "long_name": "potential temperature"}
        ds_to_expand[var] = (
            ds_to_expand["t"] * (ds_to_expand["p_f"] * p_ref_inv) ** rd_over_cp
        )
    elif var == "rho":
        attr_dict = {"units": "kg m**-3", "long_name": "density"}
        ds_to_expand[var] = ds_to_expand["p_f"] / (
            rd * ds_to_expand["t"] * (1.0 + rv_over_rd_minus_one * ds_to_expand["q"])
        )
    elif var == "w_pressure_corr":
        attr_dict = {
            "units": ds_to_expand["w"].units,
            "long_name": "Corrected pressure vertical velocity",
        }
        ds_to_expand[var] = ds_to_expand["w"] - ds_to_expand["w"][
            :, [-1], :, :
        ].values * cos_transition(
            ds_to_expand["p_f"][:, :, :, :].values,
            settings_dictionary["w_cutoff_start"],
            settings_dictionary["w_cutoff_end"],
        )
    elif var == "w_corr":
        attr_dict = {"units": "m s**-1", "long_name": "Corrected vertical velocity"}
        ds_to_expand[var] = -ds_to_expand["w_pressure_corr"] / (
            rg * ds_to_expand["rho"]
        )
    else:
        raise NotImplementedError("Variable not implemented")
    ds_to_expand[var] = ds_to_expand[var].assign_attrs(**attr_dict)


def add_auxiliary_variables(ds_to_expand, list_of_vars, settings_dictionary):
    """Wrapper for auxiliary variable calculation"""
    for var in list_of_vars:
        add_auxiliary_variable(ds_to_expand, var, settings_dictionary)


def era5_single_point(ds_domain, dictionary):
    """Extracts a local profile at the nearest point"""
    ds_at_location = ds_domain.sel(
        latitude=dictionary["lat"],
        longitude=longitude_set_meridian(dictionary["lon"]),
        method="nearest",
    )
    return ds_at_location


def era5_interp_column(ds_domain, lat_to_interp, lon_to_interp):
    """Returns the dataset interpolated to given latitude and longitude
    with latitude and longitude dimensions retained"""
    ds_at_location = ds_domain.interp(
        latitude=[lat_to_interp], longitude=[longitude_set_meridian(lon_to_interp)],
    )
    return ds_at_location


def era5_time_interp_column(ds_domain, time_to_interp, lat_to_interp, lon_to_interp):
    """Returns the dataset interpolated to given time, latitude and longitude
    with latitude and longitude dimensions retained"""
    ds_at_location = ds_domain.interp(
        time=[time_to_interp],
        latitude=[lat_to_interp],
        longitude=[longitude_set_meridian(lon_to_interp)],
    )
    return ds_at_location


def era5_weighted(ds_to_weigh, dictionary):
    """Adds weights to dictionary"""
    if "weights" in dictionary:
        if dictionary["weights"] == "area":
            ds_to_weigh.weigths = np.cos(np.deg2rad(ds_to_weigh.lat_meshgrid))
        else:
            raise Exception("weight strategy not implemented")


def era5_box_mean(ds_box, dictionary):
    """
    Calculates mean over a data_set.
    - Only use columns where the first level is higher than the local first level
    in the location of interest
    - Option to weight by box size?
    """
    era5_weighted(ds_box, dictionary)
    if "mask" in dictionary:
        mask = era5_mask(ds_box, dictionary)
        ds_mean = ds_box.where(mask).mean(("latitude", "longitude"), keep_attrs=True)
    else:
        ds_mean = ds_box.mean(("latitude", "longitude"), keep_attrs=True)
    return ds_mean


def era5_add_lat_lon_meshgrid(ds_to_extend):
    """Adds a [lat, lon] meshgrid to a dataset, useful for gradients in order to work
    around nan values on edge"""
    lon_mesh, lat_mesh = np.meshgrid(ds_to_extend.longitude, ds_to_extend.latitude)
    ds_to_extend["lon_meshgrid"] = (
        ("latitude", "longitude"),
        lon_mesh,
        {
            "long_name": "longitude on meshgrid",
            "units": ds_to_extend["longitude"].units,
        },
    )
    ds_to_extend["lat_meshgrid"] = (
        ("latitude", "longitude"),
        lat_mesh,
        {"long_name": "latitude on meshgrid", "units": ds_to_extend["latitude"].units},
    )


def calc_lat_lon_angle(lat1, lon1, lat2, lon2):
    """Calculates angle given pairs of latitude and longitude in radians"""
    dlon = lon2 - lon1
    lat_lon_angle = np.arctan2(
        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon),
        np.sin(dlon) * np.cos(lat2),
    )
    return lat_lon_angle


def era5_adv_tendencies(ds_profile, list_of_vars, dictionary):
    """Add variables defined in list to dictionary"""
    ds_out = xr.Dataset(coords={"time": ds_profile.time, "lev": ds_profile.lev})
    for variable in list_of_vars:
        tendency_array = (
            (ds_profile["u"].values - dictionary["u_traj"])
            * ds_profile["d" + variable + "dx"].values
            + (ds_profile["v"].values - dictionary["v_traj"])
            * ds_profile["d" + variable + "dy"].values
        )
        ds_out[variable + "_advtend"] = (
            ("time", "lev"),
            tendency_array,
            {
                "long_name": ds_profile[variable].long_name + " tendency (advection)",
                "units": ds_profile[variable].units + " s**-1",
            },
        )
        if dictionary["gradients_strategy"] == "both":
            tendency_array = (
                (ds_profile["u"].values - dictionary["u_traj"])
                * ds_profile["d" + variable + "dx_bound"].values
                + (ds_profile["v"].values - dictionary["v_traj"])
                * ds_profile["d" + variable + "dy_bound"].values
            )
            ds_out[variable + "_advtend_bound"] = (
                ("time", "lev"),
                tendency_array,
                {
                    "long_name": ds_profile[variable].long_name
                    + " tendency (advection, boundaries)",
                    "units": ds_profile[variable].units + " s**-1",
                },
            )

    return ds_out


def add_geowind_around_centre(ds_profile, dictionary):
    """Calculates the geostrophic wind at the centre of the box, based
    on mean profiles of density and gradients of pressure"""
    lat_centre = dictionary["lat"]
    f_cor = 2.0 * omega * np.sin(np.deg2rad(lat_centre))
    u_geo = -(1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdy"]
    v_geo = (1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdx"]
    ds_profile["ug"] = (
        ("time", "lev"),
        u_geo,
        {"long_name": "U component of geostrophic wind", "units": "m s**-1",},
    )
    ds_profile["vg"] = (
        ("time", "lev"),
        v_geo,
        {"long_name": "V component of geostrophic wind", "units": "m s**-1",},
    )
    if dictionary["gradients_strategy"] == "both":
        u_geo_bound = -(1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdy_bound"]
        v_geo_bound = (1.0 / (f_cor * ds_profile["rho"])) * ds_profile["dp_fdx_bound"]
        ds_profile["ug_bound"] = (
            ("time", "lev"),
            u_geo_bound,
            {
                "long_name": "U component of geostrophic wind (boundaries)",
                "units": "m s**-1",
            },
        )
        ds_profile["vg_bound"] = (
            ("time", "lev"),
            v_geo_bound,
            {
                "long_name": "V component of geostrophic wind (boundaries)",
                "units": "m s**-1",
            },
        )


def add_globals_attrs_to_ds(ds_to_add_to):
    """Adds global attributes to datasets"""
    global_attrs = {
        r"Conventions": r"CF-1.7",
        r"Contact": r"l.c.denby[at]leeds[dot]ac[dot again]uk s.boeing[at]leeds[dot]ac[dot again]uk",
        r"ERA5 reference": r"Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., ... & Simmons, A. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society.",
        r"Created": datetime.datetime.now().isoformat(),
        r"Created with": r"https://github.com/EUREC4A-UK/lagtraj",
    }
    for attribute in global_attrs:
        ds_to_add_to.attrs[attribute] = global_attrs[attribute]


def fix_units(ds_to_fix):
    """Changes units of ERA5 data to make them compatible with the cf-checker"""
    units_dict = {
        "(0 - 1)": "1",
        "m of water equivalent": "m",
        "~": "1",
    }
    for variable in ds_to_fix.variables:
        if hasattr(ds_to_fix[variable], "units"):
            these_units = ds_to_fix[variable].units
            if these_units in units_dict:
                ds_to_fix[variable].attrs["units"] = units_dict[these_units]


def dummy_forcings(mf_dataset, forcings_dict):
    """Forcings example: needs to be integrated into main functionality"""
    ds_out = xr.Dataset()
    ds_traj = xr.open_dataset(forcings_dict["traj_file"])
    for index in range(len(ds_traj["time"])):
        this_time = ds_traj["time"][index]
        # Ugly
        mf_index = np.argmax(mf_dataset["time"] == this_time).values
        ds_time = mf_dataset.isel(time=[mf_index])
        half_averaging_width = 0.5 * forcings_dict["averaging_width"]
        lats_lons_dict = {
            "lat_min": ds_traj["lat_traj"][index].values - half_averaging_width,
            "lat_max": ds_traj["lat_traj"][index].values + half_averaging_width,
            "lon_min": longitude_set_meridian(ds_traj["lon_traj"][index].values)
            - half_averaging_width,
            "lon_max": longitude_set_meridian(ds_traj["lon_traj"][index].values)
            + half_averaging_width,
            "lat": ds_traj["lat_traj"][index].values,
            "lon": longitude_set_meridian(ds_traj["lon_traj"][index].values),
            "u_traj": ds_traj["u_traj"][index].values,
            "v_traj": ds_traj["v_traj"][index].values,
        }
        lats_lons_dict.update(forcings_dict)
        out_levels = np.arange(0, 10000.0, 40.0)

        ds_smaller = ds_full.sel(
            latitude=slice(dictionary["lat_max"], dictionary["lat_min"]),
            longitude=slice(
                longitude_set_meridian(dictionary["lon_min"]),
                longitude_set_meridian(dictionary["lon_max"]),
            ),
        )

        add_heights_and_pressures(ds_smaller)
        add_auxiliary_variables(
            ds_smaller, ["theta", "rho", "w_pressure_corr", "w_corr"], lats_lons_dict
        )
        ds_time_height = era5_on_height_levels(ds_smaller, out_levels)
        era5_add_lat_lon_meshgrid(ds_time_height)
        ds_profiles = ds_domain.sel(
            latitude=dictionary["lat"],
            longitude=longitude_set_meridian(dictionary["lon"]),
            method="nearest",
        )

        ds_era5_mean = era5_box_mean(ds_time_height, lats_lons_dict)
        for variable in ds_era5_mean.variables:
            if variable not in ["time", "lev"]:
                ds_profiles[variable + "_mean"] = ds_era5_mean[variable]
        ds_gradients = era5_gradients(
            ds_time_height, ["u", "v", "p_f", "theta", "q"], lats_lons_dict
        )
        ds_time_step = xr.merge((ds_gradients, ds_profiles))
        ds_tendencies = era5_adv_tendencies(
            ds_time_step, ["u", "v", "theta", "q"], lats_lons_dict
        )
        ds_time_step = xr.merge((ds_time_step, ds_tendencies))
        add_geowind_around_centre(ds_time_step, lats_lons_dict)
        ds_time_step.reset_coords(["latitude", "longitude"])
        ds_out = xr.combine_by_coords((ds_out, ds_time_step))
    # Add trajectory information
    ds_out = xr.combine_by_coords((ds_out, ds_traj))
    fix_units(ds_out)
    add_globals_attrs_to_ds(ds_out)
    ds_out.attrs.update(forcings_dict)
    ds_out.to_netcdf("ds_along_traj.nc")
