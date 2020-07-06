import numpy as np
from scipy.constants import pi

from .integration.constants import r_earth
from ..utils import geometry


def extrapolate_posn_with_fixed_velocity(lat, lon, u_vel, v_vel, dt):
    """Extrapolates a position from point (lat, lon) given velocities
    (u_vel, v_vel) (in m/s) and a duration dt (in seconds)"""
    if dt < 0.0:
        raise Exception("Expecting positive dt in extrapolation")
    theta = np.arctan2(v_vel, u_vel) % (2 * pi)
    bearing = ((pi / 2) - theta) % (2 * pi)
    dist = np.sqrt(u_vel ** 2 + v_vel ** 2) * dt
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    traced_lat_rad = np.arcsin(
        np.sin(lat_rad) * np.cos(dist / r_earth)
        + np.cos(lat_rad) * np.sin(dist / r_earth) * np.cos(bearing)
    )
    traced_lon_rad = lon_rad + np.arctan2(
        np.sin(bearing) * np.sin(dist / r_earth) * np.cos(lat_rad),
        np.cos(dist / r_earth) - np.sin(lat_rad) * np.sin(lat_rad),
    )
    traced_lat = np.rad2deg(traced_lat_rad)
    traced_lon = geometry.longitude_set_meridian(np.rad2deg(traced_lon_rad))
    return traced_lat, traced_lon


def trace_two_way(lat, lon, u_vel, v_vel, dt_forward, dt_total):
    """Calculates both previous and next position given a point in between"""
    dt_backward = dt_total - dt_forward

    forward_lat, forward_lon = trace_one_way(lat, lon, u_vel, v_vel, dt_forward)
    backward_lat, backward_lon = trace_one_way(lat, lon, -u_vel, -v_vel, dt_backward)
    return forward_lat, forward_lon, backward_lat, backward_lon
