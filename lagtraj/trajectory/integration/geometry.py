


def longitude_set_meridian(longitude):
    """Sets longitude to be between -180 and 180 degrees"""
    return (longitude + 180.0) % 360.0 - 180.0
