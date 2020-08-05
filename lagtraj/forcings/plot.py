if __name__ == "__main__":  # noqa
    import matplotlib

    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import xarray as xr

from .profile_calculation import FORCING_VARS


def main(ds):
    N_vars = len(FORCING_VARS)
    fig, axes = plt.subplots(nrows=N_vars, figsize=(10, 3 * N_vars), sharex=True)

    for v, ax in zip(FORCING_VARS, axes):
        v_adv = f"d{v}dt_adv"
        ds[v_adv].plot(ax=ax, y="level")

    vel_method_details = ""
    if hasattr(ds, "velocity_method_kwargs_height"):
        vel_method_details = (
            f"velocity method using {ds.velocity_method_kwargs_height}m height"
        )

    plt.suptitle(
        f"{ds.name} {ds.trajectory_type} trajectory\n{ds.domain_name} domain\n"
        f"`{ds.velocity_method}` {vel_method_details}\n",
        y=1.01,
    )

    return fig


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("forcing_profile.nc")
    args = argparser.parse_args()

    input_filename = vars(args)["forcing_profile.nc"]
    ds = xr.open_dataset(input_filename)
    fig = main(ds=ds)

    output_filename = input_filename.replace(".nc", ".png")
    plt.savefig(output_filename, fig=fig, bbox_inches="tight")
    print(f"Saved plot to {output_filename}")
