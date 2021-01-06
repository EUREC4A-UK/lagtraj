from pathlib import Path


from .. import DEFAULT_ROOT_DATA_PATH
from ..forcings.load import load_data as load_forcing_data
from . import process
from ..utils import optional_debugging


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("forcing")
    argparser.add_argument("target_name")
    argparser.add_argument(
        "-d", "--data-path", default=DEFAULT_ROOT_DATA_PATH, type=Path
    )
    argparser.add_argument("--debug", default=False, action="store_true")
    args = argparser.parse_args()

    try:
        ds_forcing = load_forcing_data(
            root_data_path=args.data_path, forcing_name=args.forcing
        )
    except FileNotFoundError:
        raise Exception(
            f"The output file for forcing `{args.forcing}`"
            " couldn't be found. Please create the forcing by running: \n"
            f"    python -m lagtraj.forcings.create {args.forcing}\n"
            "and then run the forcing creation again"
        )

    with optional_debugging(args.debug):
        output_file_path = process.export_for_target(
            ds_forcing=ds_forcing,
            target_name=args.target_name,
            root_data_path=args.data_path,
        )

    print("Wrote converted forcing file to `{}`".format(output_file_path))


if __name__ == "__main__":
    main()
