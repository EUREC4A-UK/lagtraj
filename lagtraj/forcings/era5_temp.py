import xarray as xr


def dummy_forcings(mf_list, forcings_dict):
    """Forcings example: needs to be integrated into main functionality"""
    ds_out = xr.Dataset()
    ds_traj = xr.open_dataset(forcings_dict["traj_file"])
    for index in range(len(ds_traj["time"])):
        this_time = ds_traj["time"][index].values
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
        ds_smaller = era5_subset_by_time(mf_list, this_time, lats_lons_dict)
        add_heights_and_pressures(ds_smaller)
        add_auxiliary_variables(
            ds_smaller, ["theta", "rho", "w_pressure_corr", "w_corr"], lats_lons_dict
        )
        ds_time_height = era5_on_height_levels(ds_smaller, out_levels)
        ds_smaller.close()
        era5_add_lat_lon_meshgrid(ds_time_height)
        ds_profiles = era5_single_point(ds_time_height, lats_lons_dict)
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
        ds_time_height.close()
        ds_out = xr.combine_by_coords((ds_out, ds_time_step))
        ds_gradients.close()
        ds_era5_mean.close()
        ds_profiles.close()
        ds_tendencies.close()
        ds_time_step.close()
        ds_out.to_netcdf("ds_along_traj.nc")
    # Add trajectory information
    ds_out = xr.combine_by_coords((ds_out, ds_traj))
    fix_units(ds_out)
    ds_out["latitude"].attrs = {"long_name": "latitude", "units": "degrees_north"}
    ds_out["longitude"].attrs = {"long_name": "longitude", "units": "degrees_east"}
    add_globals_attrs_to_ds(ds_out)
    add_dict_to_global_attrs(ds_out, forcings_dict)
    ds_out.to_netcdf("ds_along_traj.nc")
