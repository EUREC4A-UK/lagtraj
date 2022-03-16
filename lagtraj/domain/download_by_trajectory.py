from .download import _run_cli


def cli(args):
    _run_cli(args=args, timedomain_lookup="by_trajectory")


if __name__ == "__main__":
    cli(args=None)
