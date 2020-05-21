

def _load_domain_data(domain):
    raise NotImplementedError


def get_available_timesteps(domain):
    return _load_domain_data(domain=domain).time


def calculate_timestep(da_pt, domain):
    """
    dt_pt.time: time for forcing
    dt_pt.lat: 
    dt_pt.lon:

    """

    domain_data = # load up data downloaded for `domain`

    if dt_pt.time is not in domain_data.time:
        raise NotImplementedError("Forcings based on era5 data cannot be"
                "interpolated between model data timesteps, and request time is"
                "outside")

    if not da_pt.levels.units == "m":
        raise InvalidLevelsDefinition(
            "ERA5 data can currently only be generated on height levels"
        )


    ds_timestep_forcing = # extract forcing profile from domain

    return ds_timestep_forcing



def append_timestep(forcings_dict, forcing_data, trajectory_data, timestep):
    if (forcings_dict["source"]).lower() == "era5":
        append_era5_timestep(forcings_dict, forcing_data, trajectory_data, timestep)


def append_era5_timestep(forcings_dict, forcing_data, trajectory_data, timestep):
    # FOR EACH TIME STEP
    # find the right input data
    # reinterpolate data to height or pressure with effective height (tbd)
    # create 'mask' for forcings on domain
    # calculate profiles and forcings
    lon_indices, lat_indices = calculate_current_lat_lon_range(
        forcings_dict, trajectory_data, timestep
    )
    timestep_input_data = xr.Dataset()
    append_single_level_an(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    append_single_level_fc(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    append_model_level_an(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    append_model_level_fc(
        timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
    )
    add_heights_and_pressures(timestep_input_data)
    timestep_interpolated_data = interpolate_to_height_levels(timestep_input_data)
    timestep_forcing_data = calculate_forcings(
        forcings_dict, timestep_interpolated_data
    )
    # append to forcing_data using xarray concatenation


def calculate_current_lat_lon_range(forcings_dict, trajectory_data, timestep):
    lon_indices = []
    lat_indices = []
    return lon_indices, lat_indices


def append_single_level_an(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def append_single_level_fc(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def append_model_level_an(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def append_model_level_fc(
    timestep_input_data, trajectory_data, timestep, lon_indices, lat_indices
):
    pass


def add_heights_and_pressures(timestep_input_data):
    pass


def interpolate_to_height_levels(timestep_input_data):
    timestep_interpolated_data = []
    return timestep_interpolated_data


def calculate_forcings(forcings_dict, forcing_interpolated_data):
    timestep_forcing_data = []
    return timestep_forcing_data


def export_to_hightune(forcing_data):
    pass


def export_to_ecmwf(forcing_data):
    pass


if __name__ == "__main__":
    main()
