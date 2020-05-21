import dateutil.parser

import yaml
import os
from pathlib import Path

from .sources import era5
from . import LatLonBoundingBox, LatLonSamplingResolution

# TODO
# - Think about a way to put multiple days in one file, which may speed up
#   downloading?
#   (on request: check for each day individually if it is there already, add
#   the missing ones to a list, download the request, split the files by date
#   again)
# - Make both NetCDF and Grib downloads an option?
# - Add a flag for putting a request online, but without download (only submit
#   a request). Then simply run again to download.


def download(dst_path, source, t_start, t_end, bbox, latlon_sampling,
             overwrite_existing=False):
    """
    Download all data from a given `source` (fx `era5`) to `dst_path` over time
    range `t_start` to `t_end`, in `bbox` with `latlon_sampling`
    """

    domain_data_path = Path(dst_path)/source.lower()

    if source.lower() == "era5":
        era5.download_data(path=domain_data_path,
                           t_start=t_start, t_end=t_end,
                           bbox=bbox, latlon_sampling=latlon_sampling,
                           overwrite_existing=overwrite_existing
                           )
    else:
        raise NotImplementedError(
            "Source type `{}` unknown. Should for example be 'era5'".format(
                source
            )
        )


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("input_file")
    argparser.add_argument("start_date", type=dateutil.parser.parse)
    argparser.add_argument("end_date", type=dateutil.parser.parse)
    argparser.add_argument(
        "-d", "--data-path", default=os.getcwd(), type=Path
    )
    argparser.add_argument(
        "-o", "--overwrite", dest="l_overwrite", action="store_true",
        default=False
    )
    args = argparser.parse_args()

    with open(args.input_file) as fh:
        domain_params = yaml.load(fh, Loader=yaml.FullLoader)

    bbox = LatLonBoundingBox(
        lat_min=domain_params['lat_min'],
        lon_min=domain_params['lon_min'],
        lat_max=domain_params['lat_max'],
        lon_max=domain_params['lon_max'],
    )

    latlon_sampling = LatLonSamplingResolution(
        lat=domain_params['lat_samp'], lon=domain_params['lon_samp']
    )

    domain_data_path = args.data_path/domain_params['name']

    download(dst_path=domain_data_path, source=domain_params['source'],
             t_start=args.start_date, t_end=args.end_date,
             bbox=bbox, latlon_sampling=latlon_sampling,
             overwrite_existing=args.l_overwrite)
