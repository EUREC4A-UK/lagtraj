import dateutil.parser

import yaml
from pathlib import Path
import sys
import textwrap

from .sources import era5
from .. import input_examples, DEFAULT_DATA_PATH
from . import LatLonBoundingBox, LatLonSamplingResolution


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
        "-d", "--data-path", default=DEFAULT_DATA_PATH, type=Path
    )
    argparser.add_argument(
        "-o", "--overwrite", dest="l_overwrite", action="store_true",
        default=False
    )
    args = argparser.parse_args()

    if args.input_file.startswith('lagtraj://'):
        try:
            input_name = (args.input_file.replace('lagtraj://', '')
                                         .replace('.yaml', ''))

            domain_params = input_examples.attempt_load(
                input_name=input_name, input_type="domain"
            )

            input_file_local = args.data_path/"domains"/input_name/"meta.yaml"
            if not input_file_local.exists():
                input_file_local.parent.mkdir(exist_ok=True, parents=True)
                with open(input_file_local, 'w') as fh:
                    fh.write(yaml.dump(domain_params))

        except input_examples.LagtrajExampleDoesNotExist:
            print("The requested domain ({}) isn't currently available"
                  "in lagtraj\n(maybe you could add it with a pull-request"
                  " at https://github.com/EUREC4A-UK/lagtraj/pulls)\n\n"
                  "The domains currently defined in lagtraj are:"
                  "".format(input_name))
            print()
            input_examples.print_available(input_types=['domains'])
            sys.exit(1)
    else:
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

    domain_data_path = args.data_path/"domains"/domain_params['name']

    download(dst_path=domain_data_path, source=domain_params['source'],
             t_start=args.start_date, t_end=args.end_date,
             bbox=bbox, latlon_sampling=latlon_sampling,
             overwrite_existing=args.l_overwrite)
