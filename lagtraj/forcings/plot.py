if __name__ == "__main__":  # noqa
    import matplotlib

    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray as xr
import numpy as np
import cartopy.crs as ccrs

from .profile_calculation import FORCING_VARS
from ..domain import load as domain_load
from ..trajectory import load as trajectory_load
from ..trajectory import plot as trajectory_plot


def main(ds):
    N_vars = len(FORCING_VARS)
    fig, axes = plt.subplots(nrows=N_vars, figsize=(10, 3 * N_vars), sharex=True)

    for v, ax in zip(FORCING_VARS, axes):
        v_adv = f"d{v}dt_adv"
        ds[v_adv].plot(ax=ax, y="level")

    title = f"{ds.name} {ds.trajectory_type} trajectory\n{ds.domain_name} domain\n"
    if hasattr(ds, "velocity_method"):
        title += (f"{ds.velocity_method} velocity method using "
                  "{ds.velocity_method_kwargs_height}m height\n")

    plt.suptitle(title, y=1.01)

    return fig


def _add_overview_axes(fig, gs):
    proj_overview = ccrs.PlateCarree()
    ax = fig.add_subplot(gs, projection=proj_overview)
    ax.coastlines(resolution="10m")
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    return ax


def plot_with_overview(
    ds,
    tn,
    forcing_vars=["dqdt_adv", "dtdt_adv"],
    domain_var="q",
    overview_window_width=4,
):
    """
    Produce a forcing plot with timestep `tn` highlighted together with
    overview plots of domain data variable `domain_var`. The width over the
    overview plot is set with `overview_window_width` (in degrees)
    """
    ds_forcing = ds
    ds_domain = domain_load.load_data(
        root_data_path="data", name=ds_forcing.domain_name
    )
    ds_traj = trajectory_load.load_data(
        root_data_path="data", name=ds_forcing.trajectory_name
    )
    N_vars = len(forcing_vars)

    figwidth = 12
    subplot_height = 3
    traj_color = "red"

    N_vars = len(forcing_vars)
    figsize = (figwidth, 4 + subplot_height * N_vars)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2 + N_vars, 2)
    ds_forc_tn = ds_forcing.isel(time=tn)

    lat0, lon0 = ds_forc_tn.lat, ds_forc_tn.lon
    domain_window = dict(
        lat=slice(lat0 - overview_window_width / 2, lat0 + overview_window_width / 2),
        lon=slice(lon0 - overview_window_width / 2, lon0 + overview_window_width / 2),
        time=ds_forc_tn.time,
    )
    da_domain = ds_domain[domain_var].sel(**domain_window).sum(dim="level")["q"]

    ax_domain = _add_overview_axes(fig=fig, gs=gs[:2, 0])
    da_domain.plot(ax=ax_domain)
    ax_satellite = _add_overview_axes(fig=fig, gs=gs[0:2, 1])
    ax_satellite.set_extent(ax_domain.get_extent())

    traj_plot_kwargs = dict(
        ds=ds_traj,
        add_ref="eurec4a_circle",
        color=traj_color,
        reference_time=ds_forc_tn.time,
    )
    trajectory_plot.main(ax=ax_domain, **traj_plot_kwargs)
    trajectory_plot.main(ax=ax_satellite, **traj_plot_kwargs)

    ax = None
    for n, v in enumerate(forcing_vars):
        ax = fig.add_subplot(gs[n + 2, :], sharex=ax)
        ds_forcing[v].plot(ax=ax, y="level")
        # .item() doesn't return a np.datetime64 object sadly, so we have to
        # make our own...
        t0 = np.datetime64(ds_forc_tn[v].time.item(), "ns")
        ax.axvline(x=t0, color="black", linestyle="--", alpha=0.5)

    fig.tight_layout()

    title = f"{ds.name} {ds.trajectory_type} trajectory\n{ds.domain_name} domain\n"
    if hasattr(ds, "velocity_method"):
        title += (
            f"{ds.velocity_method} velocity method using "
            "{ds.velocity_method_kwargs_height}m height\n"
        )

    plt.suptitle(title, y=1.01)

    return fig


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("forcing_profile.nc")
    argparser.add_argument("--highlight-timestep", type=int, nargs="*")
    args = argparser.parse_args()

    input_filename = vars(args)["forcing_profile.nc"]
    ds = xr.open_dataset(input_filename)

    if args.highlight_timestep:
        for tn in args.highlight_timestep:
            fig = plot_with_overview(ds=ds, tn=tn)
            output_filename = input_filename.replace(".nc", f"_tn{tn}.png")
            plt.savefig(output_filename, fig=fig, bbox_inches="tight")
            print(f"Saved plot to {output_filename}")
    else:
        fig = main(ds=ds)
        output_filename = input_filename.replace(".nc", ".png")
        plt.savefig(output_filename, fig=fig, bbox_inches="tight")
        print(f"Saved plot to {output_filename}")
