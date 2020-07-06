import tarfile
import numpy as np

from lagtraj.domain import load, download as domain_download
from lagtraj.trajectory import load as trajectory_load
from lagtraj import DEFAULT_ROOT_DATA_PATH


def main(
    trajectory_name="lagtraj://eurec4a_20191209_12_lin",
    output_filename="lagtraj.testdata.tar.gz",
):
    p_root = DEFAULT_ROOT_DATA_PATH

    trajectory_defn = trajectory_load.load_definition(p_root, name=trajectory_name)
    t_min = trajectory_defn.origin.datetime - trajectory_defn.duration.backward
    t_max = trajectory_defn.origin.datetime + trajectory_defn.duration.forward
    domain_name = trajectory_defn.domain

    domain_download.download_named_domain(
        data_path=DEFAULT_ROOT_DATA_PATH,
        name=domain_name,
        start_date=t_min,
        end_date=t_max,
    )
    if not domain_download.download_complete(data_path=p_root, domain_name=domain_name):
        raise Exception(
            "Some data still isn't downloaded, check the download"
            " requests have completed and run again"
        )
    ds = load.load_data(DEFAULT_ROOT_DATA_PATH, name=domain_name)

    assert not np.isnan(ds.isel(time=0).u.mean())

    p = p_root / "domains" / f"{domain_name}_data"
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(p, arcname=p.name)

    print(f"Domain data for {domain_name} written to {output_filename}")


if __name__ == "__main__":
    main()
