from collections import namedtuple


LatLonBoundingBox = namedtuple("LatLonBoundingBox",
                               ["lat_min", "lat_max", "lon_min", "lon_max"]
                               )

LatLonSamplingResolution = namedtuple("LatLonSamplingResolution",
                                      ["lat", "lon"]
                           )
