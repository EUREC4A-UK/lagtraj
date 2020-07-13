import numpy as np

from .constants import r_earth, pi


def lat_dist(lat1, lat2):
    """Distance for difference in latitude only, at constant longitude"""
    return np.deg2rad(lat2 - lat1) * r_earth


def lon_dist(lon1, lon2, lat):
    """Distance for difference in longitude, for given latitude"""
    dlon_rad = np.deg2rad(lon2 - lon1)
    dlon_rad = (dlon_rad + pi) % (2 * pi) - pi
    return np.cos(np.deg2rad(lat)) * dlon_rad * r_earth


def calc_haver_dist(lat1, lon1, lat2, lon2):
    """Calculates distance given pairs of latitude and longitude in radians"""
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    haver = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    arc_dist = 2 * np.arctan2(np.sqrt(haver), np.sqrt(1.0 - haver))
    haver_dist = r_earth * arc_dist
    return haver_dist
