import numpy as np

from .. import interpolation
from ...domain import interpolate_to_height_levels, interpolate_to_pressure_levels


def weighted_velocity(ds_column, pres_cutoff_start, pres_cutoff_end):
    """Weighted velocity: needs more work"""
    height_factor = interpolation.cos_transition(
        ds_column["p_f"][:, 1:, :, :].values, pres_cutoff_start, pres_cutoff_end
    )
    weights = (
        (ds_column["p_h"][:, :-1, :, :].values - ds_column["p_h"][:, 1:, :, :].values)
        * ds_column["q"][:, 1:, :, :].values
        * height_factor
    )
    inv_weights = 1.0 / np.sum(weights)
    u_weighted = inv_weights * np.sum(ds_column["u"][:, 1:, :, :].values * weights)
    v_weighted = inv_weights * np.sum(ds_column["v"][:, 1:, :, :].values * weights)
    return u_weighted, v_weighted


def velocity_at_height(ds_column, height):
    """Velocit at one height: needs more work"""
    ds_on_height_level = interpolate_to_height_levels(ds_column, height)
    # For a single colum, data dimensions are all 1
    return np.mean(ds_on_height_level["u"]), np.mean(ds_on_height_level["v"])


def velocity_at_pressure(ds_column, pressure):
    """Velocit at one height: needs more work"""
    ds_on_pressure_level = interpolate_to_pressure_levels(ds_column, pressure)
    # For a single colum, data dimensions are all 1
    return np.mean(ds_on_pressure_level["u"]), np.mean(ds_on_pressure_level["v"])


def estimate_horizontal_velocities(ds_column, method, **kwargs):
    """Estimate the zonal and meridonal winds using a specific method"""
    if method == "lower_troposphere_humidity_weighted":
        if "pres_cutoff_start" not in kwargs or "pres_cutoff_end" not in kwargs:
            raise Exception(
                f"To use the `{method}` velocity method the"
                " `pres_cutoff_start` and `pres_cutoff_end` kwargs"
                " are required"
            )

        u_traj, v_traj = weighted_velocity(ds_column, **kwargs)
    elif method == "single_height_level":
        if "height" not in kwargs:
            raise Exception(
                f"To use the `{method}` velocity method the"
                " `height` kwarg is required"
            )
        u_traj, v_traj = velocity_at_height(ds_column, **kwargs)
    elif method == "single_pressure_level":
        if "pressure" not in kwargs:
            raise Exception(
                f"To use the `{method}` velocity method the"
                " `pressure` kwarg is required"
            )
        u_traj, v_traj = velocity_at_pressure(ds_column, **kwargs)
    elif method == "column_mean":
        u_traj = np.float64(ds_column.u.mean())
        v_traj = np.float64(ds_column.v.mean())
    else:
        raise NotImplementedError(
            f"`{method}` trajectory velocity method" " not implemented"
        )

    return np.float64(u_traj), np.float64(v_traj)
