from collections import namedtuple
import numpy as np


from ..utils.gradient_calculation import calc_horizontal_gradients
from ..domain import calc_auxiliary_variable as calc_auxiliary_domain_variable
from ..domain import interpolate_to_height_levels as interpolate_domain_to_height_levels
from ..domain.mask import calc_mask


# list of scalars we want to compute forcings of, TODO: move into yaml input
# definitions
FORCING_VARS = [
    "u",
    "v",
    "theta",
    "q",
    "clwc",
    "ciwc",
    "t",
    "p_f",
    "t_l",
    "cc",
    "q_t",
    "r_t",
    "r_v",
    "r_l",
    "r_i",
    "theta_l",
]

# list of variables which must be derived, TODO: move itno yaml input
# definition, NOTE: the order is important here, `rho` must be derived first as
# it's needed for `w_corr` and `w_pressure_corr`, TODO: the order here should
# be automatically resolved in a module specialied to era5 (so we can
# generalise to model input from other models)
AUXILIARY_VARS = [
    "theta",
    "rho",
    "w_pressure_corr",
    "w_corr",
    "t_l",
    "q_t",
    "q_t_hydromet",
    "r_t",
    "r_v",
    "r_l",
    "r_i",
    "theta_l",
]

# Geostrophic winds can only be calculated after the profile calculations
FINAL_VARS = ["u_g", "v_g"]

ForcingSamplingDefinition = namedtuple(
    "ForcingSamplingDefinition",
    ["gradient_method", "averaging_width", "time_sampling_method", "mask"],
)


def _build_domain_profile(da_field, method="single_point", **kwargs):
    """
    Construct from `da_field` a representative vertical profile using `method`
    """
    if method == "single_point":
        if "reference_point" not in kwargs:
            raise Exception(
                "using method {method} profile method requires the"
                " `reference_point` kwarg to be provided"
            )
        # we .squeeze() here to ensure we don't try interpolate on dimension
        # which aren't necessary
        return da_field.squeeze().interp(**kwargs["reference_point"].squeeze())
    elif method == "mean":
        # have to provide `dtype` kwarg otherwise `bottleneck` might use
        # float32 to calculate means
        return da_field.mean(dim=("lat", "lon"), dtype=np.float64)
    else:
        raise NotImplementedError(method)


def compute_adv_profile(ds_profile, da_domain, gradient_method):
    """
    Compute the horizontal advective tendency profile and reference profile for
    a field `da_domain` (named `phi` here) for a point on trajectory (defined
    by `lat` and `lon` on `ds_profile`) with horizontal velocity `ds_traj.u`
    and `ds_traj.v`

        dphi/dt = - dphi/dx * u - dphi/dy * v

    given the vertical `ds_profile` containing ambient horizontal wind
    components (u, v) and horizontal gradients (dphi/dx, dphi/dy)
    """
    # build a reference point for the interpolations required when calculating
    # the representative vertical profile and horizontal gradients
    ds_ref_pt = ds_profile[["lat", "lon", "time"]]

    da_dphidx, da_dphidy = calc_horizontal_gradients(
        da_field=da_domain, method=gradient_method, ds_ref_pt=ds_ref_pt
    )

    # compute the relative velocities (`_local` being the velocities
    # interpolated to the trajectory (lat, lon) and `_traj` being the velocity
    # of the trajectory itself)
    da_u_rel = ds_profile.u_local - ds_profile.u_traj
    da_v_rel = ds_profile.v_local - ds_profile.v_traj

    # dphi/dt = - dphi/dx * u - dphi/dy * v
    da_dphidt = -da_dphidx * da_u_rel - da_dphidy * da_v_rel

    da_dphidt.attrs["long_name"] = f"{da_domain.long_name} tendency (advection)"
    da_dphidt.attrs["units"] = f"{da_domain.units} s**-1"

    da_dphidt = da_dphidt.squeeze()

    return da_dphidt, da_dphidx, da_dphidy


class InvalidLevelsDefinition(Exception):
    pass


def _reset_lat_lon(ds):
    # Function to ensure lat and lon are not a dimension on this data.
    ds=ds.reset_coords(["lat", "lon"])
    ds["lat"].attrs = {"long_name": "latitude", "units": "degrees_north"}
    ds["lon"].attrs = {"long_name": "longitude", "units": "degrees_east"}
    return ds


def _construct_subdomain(
    ds_profile_posn,
    ds_domain,
    latlon_sampling_window,
    mask_type,
    use_hacky_variable_preselection=False,
):
    l_lat = l_lon = latlon_sampling_window
    sampling_window = dict(
        lat=slice(ds_profile_posn.lat - l_lat / 2.0, ds_profile_posn.lat + l_lat / 2.0),
        lon=slice(ds_profile_posn.lon - l_lon / 2.0, ds_profile_posn.lon + l_lon / 2.0),
        time=ds_profile_posn.time,
    )

    # TODO: make it possible to select which variables to include in the output
    # (mean and local) profiles, for now use all variables available in the
    # domain dataset
    required_vars = list(ds_domain.data_vars)

    # clip the domain to the sampling window
    ds_subdomain = ds_domain[required_vars].sel(**sampling_window)

    # and apply a mask (we might for example only want to consider columns
    # which are over the ocean)
    da_mask = calc_mask(ds=ds_subdomain, mask_type=mask_type)
    ds_subdomain["mask"] = da_mask
    ds_subdomain = ds_subdomain.where(ds_subdomain.mask, other=np.nan)
    ds_subdomain.attrs["data_source"] = ds_domain.attrs.get("data_source")

    # TODO: this interpolates all domain variables to height levels, but we
    # should only interpolate the variables necessary to produce the requested
    # output variables (and their forcings)
    ds_subdomain_hl = interpolate_domain_to_height_levels(
        ds=ds_subdomain, height=ds_profile_posn.level
    )

    # attempt to calculate the missing variables as auxiliary variables
    for v in AUXILIARY_VARS:
        aux_kwargs = {}
        # TODO: these parameters should got into the `sampling_method` definition
        if v == "w_pressure_corr":
            aux_kwargs["w_cutoff_start"] = 70000.0
            aux_kwargs["w_cutoff_end"] = 40000.0
        ds_subdomain_hl[v] = calc_auxiliary_domain_variable(
            ds=ds_subdomain_hl, v=v, **aux_kwargs
        )

    # ensure all data is loaded into memory
    return ds_subdomain_hl.compute()


def calculate_timestep(ds_profile_posn, ds_domain, sampling_method):
    """
    given the domain data in `ds_domain` calculate the forcing profile at
    points in `ds_profile_posn`, defining the `time`, `lat` and `lon`
    and `level` (height) positions.
    """
    ds_domain_timestep = ds_domain.time.sel(time=ds_profile_posn.time)
    if ds_domain_timestep.time.count() != 1:
        raise NotImplementedError(
            "Forcings based on era5 data cannot be"
            "interpolated between model data timesteps, and requested time is"
            "outside"
        )

    if not ds_profile_posn.level.units == "m":
        raise InvalidLevelsDefinition(
            "ERA5 data can currently only be generated on height levels"
        )

    # remove the `u_traj` and `v_traj` from `ds_profile_posn` before we start
    # doing any interpolation
    u_traj, v_traj = ds_profile_posn.u_traj, ds_profile_posn.v_traj
    ds_profile_posn = ds_profile_posn.drop(["u_traj", "v_traj"])

    # extract a sampling domain ensuring that the required variables are
    # present and on the correct grid
    ds_subdomain = _construct_subdomain(
        ds_profile_posn=ds_profile_posn,
        ds_domain=ds_domain,
        latlon_sampling_window=sampling_method.averaging_width,
        mask_type=sampling_method.mask,
    )

    # start with a profile with just the horizontal wind profiles estimated at
    # the trajectory point (these are needed to compute the advective derivatives)
    ds_profile = ds_profile_posn.copy().set_coords(["time", "level", "lat", "lon"])
    ds_profile.attrs["data_source"] = "era5"
    ds_profile["u_traj"] = u_traj
    ds_profile["v_traj"] = v_traj

    # compute mean profiles and profile at trajectory (lat, lon)-point of all
    # 3D variables. For this we need a reference point for the interpolations
    # required when calculating the representative vertical profile and
    # horizontal gradients
    ds_ref_pt = ds_profile_posn[["lat", "lon", "time"]]
    for v in ds_subdomain.data_vars:
        da_field = ds_subdomain[v]
        # have to provide `dtype` kwarg otherwise `bottleneck` might use
        # float32 to calculate means
        da_v__mean = da_field.mean(
            dim=("lat", "lon"), dtype=np.float64, keep_attrs=True
        )
        da_v__mean.attrs["long_name"] = f"sampling-domain mean {da_field.long_name}"
        ds_profile[f"{v}_mean"] = da_v__mean

        da_v__local = da_field.squeeze().interp(ds_ref_pt)
        da_v__local.attrs["long_name"] = f"trajectory-centered {da_field.long_name}"
        da_v__local.attrs["units"] = da_field.units
        ds_profile[f"{v}_local"] = da_v__local

    for v in FORCING_VARS:
        da_subdomain = ds_subdomain[v]
        da_adv_profile, da_dvdx, da_dvdy = compute_adv_profile(
            da_domain=da_subdomain,
            ds_profile=ds_profile,
            gradient_method=sampling_method.gradient_method,
        )
        ds_profile[f"d{v}dt_adv"] = da_adv_profile
        ds_profile[f"d{v}dx"] = da_dvdx
        ds_profile[f"d{v}dy"] = da_dvdy

    for v in FINAL_VARS:
        aux_kwargs = {}
        ds_profile[v] = calc_auxiliary_domain_variable(ds=ds_profile, v=v, **aux_kwargs)

    ds_profile = _reset_lat_lon(ds_profile)
    return ds_profile
