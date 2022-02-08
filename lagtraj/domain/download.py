import datetime
from pathlib import Path
from time import sleep

import dateutil.parser

from .. import DEFAULT_ROOT_DATA_PATH
from ..trajectory.load import load_definition as load_traj_definition
from . import LatLonBoundingBox, LatLonSamplingResolution, build_domain_data_path
from .load import load_definition
from .sources import era5


def download(
    data_path,
    source,
    start_date,
    end_date,
    bbox,
    latlon_sampling,
    version,
    overwrite_existing=False,
):
    """
    Download all data from a given `source` (fx `era5`) to `data_path` over time
    range `start_date` to `end_date`, in `bbox` with `latlon_sampling` and store it
    with `version` as the version of this domain.
    """

    if source.lower() == "era5":
        era5.download_data(
            path=data_path,
            start_date=start_date,
            end_date=end_date,
            bbox=bbox,
            latlon_sampling=latlon_sampling,
            overwrite_existing=overwrite_existing,
            version=version,
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
        start_date=start_date,
        end_date=end_date,
        bbox=bbox,
        latlon_sampling=latlon_sampling,
        overwrite_existing=overwrite_existing,
        version=domain_params["version"],
    )


def _find_missing_files(root_data_path, domain_name, start_date, end_date):
    """
    For the selected data source return a list of the files that haven't yet been downloaded
    """
    domain_params = load_definition(domain_name=domain_name, data_path=root_data_path)
    source = domain_params["source"]

    domain_data_path = build_domain_data_path(
        root_data_path=root_data_path, domain_name=domain_params["name"]
    )

    bbox = LatLonBoundingBox(
        lat_min=domain_params["lat_min"],
        lon_min=domain_params["lon_min"],
        lat_max=domain_params["lat_max"],
        lon_max=domain_params["lon_max"],
    )

    latlon_sampling = LatLonSamplingResolution(
        lat=domain_params["lat_samp"], lon=domain_params["lon_samp"]
    )

    if source.lower() == "era5":
        return era5.find_missing_files(
            path=domain_data_path,
            start_date=start_date,
            end_date=end_date,
            bbox=bbox,
            latlon_sampling=latlon_sampling,
        )
    else:
        raise NotImplementedError(
            "Source type `{}` unknown. Should for example be 'era5'".format(source)
        )


def download_complete(root_data_path, domain_name, start_date, end_date):
    """
    Check that all files have been downloaded and that they contain the data in
    the expected date range
    """
    missing_files = _find_missing_files(
        root_data_path=root_data_path,
        domain_name=domain_name,
        start_date=start_date,
        end_date=end_date,
    )
    return len(missing_files) == 0


def list_files_still_to_download(root_data_path, domain_name, start_date, end_date):
    """
    Print out the files that haven't yet been downloaded
    """
    missing_files = _find_missing_files(
        root_data_path=root_data_path,
        domain_name=domain_name,
        start_date=start_date,
        end_date=end_date,
    )
    if len(missing_files) == 0:
        print("All files downloaded!")
    else:
        print("The following datafiles would be requested for download:")
        for item in missing_files:
            print(f"\t{item['source']}: {item['filepath'].name}")


def _parse_date(s):
    return dateutil.parser.parse(s).date()


def _make_cli_argparser(timedomain_lookup):
    import argparse

    argparser = argparse.ArgumentParser()
    if timedomain_lookup == "by_arguments":
        argparser.add_argument("domain")
        argparser.add_argument("start_date", type=_parse_date)
        argparser.add_argument("end_date", type=_parse_date)
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
    argparser.add_argument(
        "--dry-run",
        default=False,
        action="store_true",
        help="Don't actually make any data requests, only list the ones that would be made",
    )

    return argparser


def _run_cli(args=None, timedomain_lookup="by_arguments"):
    argparser = _make_cli_argparser(timedomain_lookup=timedomain_lookup)
    args = argparser.parse_args(args=args)

    if timedomain_lookup == "by_arguments":
        domain = args.domain
        start_date = args.start_date
        end_date = args.end_date
    elif timedomain_lookup == "by_trajectory":
        traj_defn = load_traj_definition(
            root_data_path=args.data_path, name=args.trajectory
        )
        domain = traj_defn.domain
        t_min = traj_defn.origin.datetime - traj_defn.duration.backward
        t_max = traj_defn.origin.datetime + traj_defn.duration.forward
        # downloads are done by date, so we just extract the date here to
        # ensure the whole day is made available
        start_date = t_min.date()
        end_date = t_max.date()
    else:
        raise NotImplementedError(timedomain_lookup)

    if args.dry_run:
        list_files_still_to_download(
            root_data_path=args.data_path,
            domain_name=domain,
            start_date=start_date,
            end_date=end_date,
        )
        return

    def _download_complete():
        return download_complete(
            args.data_path,
            domain_name=domain,
            start_date=start_date,
            end_date=end_date,
        )

    def attempt_download():
        download_named_domain(
            data_path=args.data_path,
            name=domain,
            start_date=start_date,
            end_date=end_date,
            overwrite_existing=args.l_overwrite,
        )

    while not _download_complete():
        attempt_download()
        if args.retry_rate is None:
            break
        else:
            t_now = datetime.datetime.now()
            t_now_s = t_now.strftime("%Y%m%dT%H%M")
            print(f"{t_now_s}: Sleeping {args.retry_rate}min...")
            sleep(args.retry_rate * 60.0)
            print("Retrying download")


def cli(args):
    _run_cli(args=args, timedomain_lookup="by_arguments")


def has_data_for_cli_command(args):
    argparser = _make_cli_argparser(timedomain_lookup="by_arguments")
    args = argparser.parse_args(args=args)

    missing_files = _find_missing_files(
        root_data_path=args.data_path,
        domain_name=args.domain,
        start_date=args.start_date,
        end_date=args.end_date,
    )
    return len(missing_files) == 0


def cli(args):
    _run_cli(args=args, timedomain_lookup="by_arguments")


def has_data_for_cli_command(args):
    argparser = _make_cli_argparser(timedomain_lookup="by_arguments")
    args = argparser.parse_args(args=args)

    missing_files = _find_missing_files(
        root_data_path=args.data_path,
        domain_name=args.domain,
        start_date=args.start_date,
        end_date=args.end_date,
    )
    return len(missing_files) == 0


if __name__ == "__main__":
    cli(args=None)
