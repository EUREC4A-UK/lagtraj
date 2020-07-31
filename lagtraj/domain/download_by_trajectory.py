from pathlib import Path


from .. import DEFAULT_ROOT_DATA_PATH
from .download import download_named_domain
from ..trajectory.load import load_definition as load_traj_definition


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("trajectory")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument(
        "-o", "--overwrite", dest="l_overwrite", action="store_true", default=False
    )
    args = argparser.parse_args()

    traj_defn = load_traj_definition(
        root_data_path=args.data_path, name=args.trajectory
    )

    t_min = traj_defn.origin.datetime - traj_defn.duration.backward
    t_max = traj_defn.origin.datetime + traj_defn.duration.forward

    download_named_domain(
        data_path=args.data_path,
        name=traj_defn.domain,
        start_date=t_min,
        end_date=t_max,
        overwrite_existing=args.l_overwrite,
    )
