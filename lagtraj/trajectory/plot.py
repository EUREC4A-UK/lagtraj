if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")  # noqa

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr

try:
    import twinotter.external.eurec4a

    HAS_TWINOTTER = True
except ImportError:
    HAS_TWINOTTER = False


def main(ds, add_ref=None, ax=None, reference_time="origin", **kwargs):
    new_axes = False
    if ax is None:
        map_proj = ccrs.Mercator()
        ax = plt.axes(projection=map_proj)
        ax.gridlines(draw_labels=True)
        ax.coastlines(resolution="10m")
        new_axes = True

    kwargs.setdefault("color", "darkgreen")
    kwargs.setdefault("marker", "o")
    kwargs["transform"] = ccrs.PlateCarree()

    if reference_time == "origin":
        t_ref = ds.origin_datetime
        lat_ref = ds.origin_lat
        lon_ref = ds.origin_lon
    else:
        t_ref = reference_time
        lat_ref = ds.sel(time=t_ref).lat
        lon_ref = ds.sel(time=t_ref).lon

    ds_backward = ds.sel(time=slice(None, t_ref))
    (line,) = ax.plot(
        ds_backward.lon,
        ds_backward.lat,
        markersize=1,
        alpha=0.4,
        **kwargs,
    )

    ds_forward = ds.sel(time=slice(t_ref, None))
    ax.plot(ds_forward.lon, ds_forward.lat, markersize=1, **kwargs)

    ax.plot(
        lon_ref,
        lat_ref,
        markersize=4,
        linestyle="",
        **kwargs,
    )

    if add_ref == "eurec4a_circle":
        if not HAS_TWINOTTER:
            raise Exception(
                "Please install `twinotter` to enable adding the" " EUREC4A circle"
            )
        twinotter.external.eurec4a.add_halo_circle(
            ax=ax, color=kwargs["color"], alpha=0.5
        )
        if new_axes:
            ax.set_extent([-64, -53, 10, 16], crs=ccrs.PlateCarree())

    DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S"
    t_start_str = ds.time.min().dt.strftime(DATETIME_FORMAT).item()
    t_end_str = ds.time.max().dt.strftime(DATETIME_FORMAT).item()
    plt.title(
        f"{ds.name} {ds.trajectory_type} trajectory\n{ds.domain_name} domain\n"
        f"from {t_start_str} to {t_end_str}"
        "\n\n"
    )

    return ax


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("trajectory.nc")
    argparser.add_argument("--add-ref", default=None)
    args = argparser.parse_args()

    input_filename = vars(args)["trajectory.nc"]
    ds = xr.open_dataset(input_filename)
    ax = main(ds=ds, add_ref=args.add_ref)

    output_filename = input_filename.replace(".nc", ".png")
    plt.savefig(output_filename)
    print(f"Saved plot to {output_filename}")
