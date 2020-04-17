

def validate_trajectory(ds_traj):
    required_fields = ["lat", "lon", "time"]
    missing_fields = list(filter(lambda f: f not in ds_traj,
                                 required_fields))

    if len(missing_fields) > 0:
        raise Exception("The provided trajectory is missing the following"
                        " fields: {}".format(", ".join(missing_fields)))


def validate_forcing_profiles(ds_forcing_profiles):
    required_fields = [
        "lat", "lon", "time", "level", "ddt__qv"
    ]
    missing_fields = list(filter(lambda f: f not in ds_forcing_profiles,
                                 required_fields))

    if len(missing_fields) > 0:
        raise Exception("The provided forcing profiles are missing the"
                        " following fields: {}".format(
                        ", ".join(missing_fields)))
