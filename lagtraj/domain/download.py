import dateutil.parser

from pathlib import Path

from .sources import era5
from .. import DEFAULT_ROOT_DATA_PATH
from . import LatLonBoundingBox, LatLonSamplingResolution, build_domain_data_path
from .load import load_definition


def download(
    dst_path, source, t_start, t_end, bbox, latlon_sampling, overwrite_existing=False
):
    """
    Download all data from a given `source` (fx `era5`) to `dst_path` over time
    range `t_start` to `t_end`, in `bbox` with `latlon_sampling`
    """

    domain_data_path = Path(dst_path) / source.lower()

    if source.lower() == "era5":
        era5.download_data(
            path=domain_data_path,
            t_start=t_start,
            t_end=t_end,
            bbox=bbox,
            latlon_sampling=latlon_sampling,
            overwrite_existing=overwrite_existing,
        )
    else:
        raise NotImplementedError(
            "Source type `{}` unknown. Should for example be 'era5'".format(source)
        )


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("domain")
    argparser.add_argument("start_date", type=dateutil.parser.parse)
    argparser.add_argument("end_date", type=dateutil.parser.parse)
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument(
        "-o", "--overwrite", dest="l_overwrite", action="store_true", default=False
    )
    args = argparser.parse_args()

    domain_params = load_definition(domain_name=args.domain, data_path=args.data_path)

    bbox = LatLonBoundingBox(
        lat_min=domain_params["lat_min"],
        lon_min=domain_params["lon_min"],
        lat_max=domain_params["lat_max"],
        lon_max=domain_params["lon_max"],
    )

    latlon_sampling = LatLonSamplingResolution(
        lat=domain_params["lat_samp"], lon=domain_params["lon_samp"]
    )

    domain_data_path = build_domain_data_path(
        root_data_path=args.data_path, domain_name=domain_params["name"]
    )

    download(
        dst_path=domain_data_path,
        source=domain_params["source"],
        t_start=args.start_date,
        t_end=args.end_date,
        bbox=bbox,
        latlon_sampling=latlon_sampling,
        overwrite_existing=args.l_overwrite,
    )
