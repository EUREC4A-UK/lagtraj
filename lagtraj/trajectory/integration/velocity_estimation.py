import numpy as np

from .. import interpolation
from .constants import rv_over_rd_minus_one


def weighted_velocity(ds_for_vel, trajectory_dict):
    """Weighted velocity: needs more work"""
    pres_cutoff_start = trajectory_dict["pres_cutoff_start"]
    pres_cutoff_end = trajectory_dict["pres_cutoff_end"]
    height_factor = interpolation.cos_transition(
        ds_for_vel["p_f"][:, 1:, :, :].values, pres_cutoff_start, pres_cutoff_end
    )
    weights = (
        (ds_for_vel["p_h"][:, :-1, :, :].values - ds_for_vel["p_h"][:, 1:, :, :].values)
        * ds_for_vel["q"][:, 1:, :, :].values
        * height_factor
    )
    inv_weights = 1.0 / np.sum(weights)
    u_weighted = inv_weights * np.sum(ds_for_vel["u"][:, 1:, :, :].values * weights)
    v_weighted = inv_weights * np.sum(ds_for_vel["v"][:, 1:, :, :].values * weights)
    return u_weighted, v_weighted


def velocity_at_height(ds_for_vel, trajectory_dict):
    """Velocit at one height: needs more work"""
    single_height_level = np.array([trajectory_dict["velocity_height"]])
    ds_on_height_level = era5_on_height_levels(ds_for_vel, single_height_level)
    # For a single colum, data dimensions are all 1
    return np.mean(ds_on_height_level["u"]), np.mean(ds_on_height_level["v"])


def get_velocity_from_strategy(ds_column, trajectory_dict):
    """wrapper routine, determine velocity according to strategy"""
    if trajectory_dict["velocity_strategy"] == "lower_troposphere_humidity_weighted":
        u_traj, v_traj = weighted_velocity(ds_column, trajectory_dict)
    elif trajectory_dict["velocity_strategy"] == "velocity_at_height":
        u_traj, v_traj = velocity_at_height(ds_column, trajectory_dict)
    else:
        raise NotImplementedError("Trajectory velocity strategy not implemented")
    return u_traj, v_traj
