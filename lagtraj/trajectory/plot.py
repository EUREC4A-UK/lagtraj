import matplotlib

matplotlib.use("Agg")  # noqa

import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

try:
    import twinotter.external.eurec4a

    HAS_TWINOTTER = True
except ImportError:
    HAS_TWINOTTER = False


def main(input_filename, add_ref=None):
    ds = xr.open_dataset(input_filename)

    map_proj = ccrs.Mercator()

    ax = plt.axes(projection=map_proj)

    plt_kwargs = dict(transform=ccrs.PlateCarree(), marker="o", color="darkgreen")

    ds_backward = ds.sel(time=slice(None, ds.origin_datetime))
    (line,) = ax.plot(
        ds_backward.lon, ds_backward.lat, markersize=1, alpha=0.4, **plt_kwargs,
    )

    ds_forward = ds.sel(time=slice(ds.origin_datetime, None))
    ax.plot(ds_forward.lon, ds_forward.lat, markersize=1, **plt_kwargs)

    ax.plot(
        ds.origin_lon, ds.origin_lat, markersize=4, linestyle="", **plt_kwargs,
    )

    ax.gridlines(draw_labels=True)
    ax.coastlines(resolution="10m")

    if add_ref == "eurec4a_circle":
        if not HAS_TWINOTTER:
            raise Exception(
                "Please install `twinotter` to enable adding the" " EUREC4A circle"
            )
        twinotter.external.eurec4a.add_halo_circle(ax=ax)
        ax.set_extent([-64, -53, 10, 16], crs=ccrs.PlateCarree())

    DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S"
    t_start_str = ds.time.min().dt.strftime(DATETIME_FORMAT).item()
    t_end_str = ds.time.max().dt.strftime(DATETIME_FORMAT).item()
    plt.title(
        f"{ds.name} {ds.trajectory_type} trajectory\n{ds.domain_name} domain\n"
        f"from {t_start_str} to {t_end_str}"
        "\n\n"
    )

    output_filename = input_filename.replace(".nc", ".png")
    plt.savefig(output_filename)
    print(f"Saved plot to {output_filename}")


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("trajectory.nc")
    argparser.add_argument("--add-ref", default=None)
    args = argparser.parse_args()

    main(input_filename=vars(args)["trajectory.nc"], add_ref=args.add_ref)
