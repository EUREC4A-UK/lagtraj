import numpy as np
from .constants import r_earth
from scipy.constants import pi


def lat_dist(lat1, lat2):
    """Distance for difference in latitude only, at constant longitude"""
    return np.deg2rad(lat2 - lat1) * r_earth


def lon_dist(lon1, lon2, lat):
    """Distance for difference in longitude, for given latitude"""
    dlon_rad = np.deg2rad(lon2 - lon1)
    dlon_rad = (dlon_rad + pi) % (2 * pi) - pi
    return np.cos(np.deg2rad(lat)) * dlon_rad * r_earth


def longitude_set_meridian(longitude):
    """Sets longitude to be between -180 and 180 degrees"""
    return (longitude + 180.0) % 360.0 - 180.0
