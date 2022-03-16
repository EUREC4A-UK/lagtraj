import tarfile

from lagtraj import DEFAULT_ROOT_DATA_PATH
from lagtraj.domain import download as domain_download
from lagtraj.domain import load
from lagtraj.forcings import load as forcing_load
from lagtraj.trajectory import load as trajectory_load
from lagtraj.utils import optional_debugging

TEST_FORCING_NAME = "lagtraj://eurec4a_20200202_first_short"


def main(
    forcing_name=TEST_FORCING_NAME,
    output_filename="lagtraj.testdata.tar.gz",
):
    p_root = DEFAULT_ROOT_DATA_PATH

    forcing_defn = forcing_load.load_definition(p_root, forcing_name=forcing_name)
    trajectory_name = forcing_defn.name

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
    if not domain_download.download_complete(
        root_data_path=p_root, domain_name=domain_name
    ):
        raise Exception(
            "Some data still isn't downloaded, check the download"
            " requests have completed and run again"
        )
    # attempt to load the data
    _ = load.load_data(DEFAULT_ROOT_DATA_PATH, name=domain_name)

    # NB: the lines below add all domain-data files in the domain data
    # directory to the tar-ball and so may include dates that aren't needed if
    # more data has been downloaded manually
    domain_name_local = domain_name.replace("lagtraj://", "")
    with tarfile.open(output_filename, "w:gz") as tar:
        p = p_root / "domains" / f"{domain_name_local}.yaml"
        tar.add(p, arcname=p.relative_to(p_root))
        p = p_root / "domains" / f"{domain_name_local}_data"
        tar.add(p, arcname=p.relative_to(p_root))

    print(f"Domain data for {domain_name} written to {output_filename}")


if __name__ == "__main__":
    with optional_debugging(True):
        main()
