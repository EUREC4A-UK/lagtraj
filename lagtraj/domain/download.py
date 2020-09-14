import dateutil.parser
from pathlib import Path
from time import sleep

from .sources import era5
from .. import DEFAULT_ROOT_DATA_PATH
from . import LatLonBoundingBox, LatLonSamplingResolution, build_domain_data_path
from .load import load_definition
from ..trajectory.load import load_definition as load_traj_definition


def download(
    data_path, source, t_start, t_end, bbox, latlon_sampling, overwrite_existing=False
):
    """
    Download all data from a given `source` (fx `era5`) to `data_path` over time
    range `t_start` to `t_end`, in `bbox` with `latlon_sampling`
    """

    if source.lower() == "era5":
        era5.download_data(
            path=data_path,
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


def download_named_domain(
    data_path, name, start_date, end_date, overwrite_existing=False
):
    domain_params = load_definition(domain_name=name, data_path=data_path)

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
        root_data_path=data_path, domain_name=domain_params["name"]
    )

    download(
        data_path=domain_data_path,
        source=domain_params["source"],
        t_start=start_date,
        t_end=end_date,
        bbox=bbox,
        latlon_sampling=latlon_sampling,
        overwrite_existing=overwrite_existing,
    )


def download_complete(root_data_path, domain_name):
    domain_params = load_definition(domain_name=domain_name, data_path=root_data_path)
    source = domain_params["source"]

    domain_data_path = build_domain_data_path(
        root_data_path=root_data_path, domain_name=domain_params["name"]
    )

    if source.lower() == "era5":
        return era5.all_data_is_downloaded(path=domain_data_path)
    else:
        raise NotImplementedError(
            "Source type `{}` unknown. Should for example be 'era5'".format(source)
        )


def _run_cli(timedomain_lookup="by_arguments"):
    import argparse

    argparser = argparse.ArgumentParser()
    if timedomain_lookup == "by_arguments":
        argparser.add_argument("domain")
        argparser.add_argument("start_date", type=dateutil.parser.parse)
        argparser.add_argument("end_date", type=dateutil.parser.parse)
    elif timedomain_lookup == "by_trajectory":
        argparser.add_argument("trajectory")
    else:
        raise NotImplementedError(timedomain_lookup)

    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument(
        "-o", "--overwrite", dest="l_overwrite", action="store_true", default=False
    )
    argparser.add_argument(
        "--retry-rate",
        default=None,
        type=float,
        help="Retry time delay (in minutes) when some files are still processing",
    )
    args = argparser.parse_args()

    if timedomain_lookup == "by_arguments":
        domain = args.domain
        t_min = args.start_date
        t_max = args.end_date
    elif timedomain_lookup == "by_trajectory":
        traj_defn = load_traj_definition(
            root_data_path=args.data_path, name=args.trajectory
        )
        domain = traj_defn.domain
        t_min = traj_defn.origin.datetime - traj_defn.duration.backward
        t_max = traj_defn.origin.datetime + traj_defn.duration.forward
    else:
        raise NotImplementedError(timedomain_lookup)

    def attempt_download():
        download_named_domain(
            data_path=args.data_path,
            name=domain,
            start_date=t_min,
            end_date=t_max,
            overwrite_existing=args.l_overwrite,
        )

    if args.retry_rate is not None:
        while True:
            attempt_download()
            if download_complete(args.data_path, domain_name=domain):
                break
            else:
                print(f"Sleeping {args.retry_rate}min...")
                sleep(args.retry_rate * 60.0)
                print("Retrying download")
    else:
        attempt_download()


if __name__ == "__main__":
    _run_cli(timedomain_lookup="by_arguments")
