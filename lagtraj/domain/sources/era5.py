"""
Routines for downloading era5 data for a given domain and storing it locally
"""
import datetime
import warnings
from pathlib import Path
import xarray as xr
import numpy as np

import netCDF4
import yaml

from ... import utils
from .cdsapi_request import RequestFetchCDSClient


DATE_FORMAT = "%Y-%m-%d"
TIME_FORMAT = "%h:%M:%s"
REPOSITORY_NAME = "reanalysis-era5-complete"
FILENAME_FORMAT = "{model_run_type}_{level_type}_{date}.nc"
DATA_REQUESTS_FILENAME = "data_requests.yaml"


def download_data(
    path, t_start, t_end, bbox, latlon_sampling, overwrite_existing=False
):

    dl_queries = _build_queries(
        t_start=t_start, t_end=t_end, bbox=bbox, latlon_sampling=latlon_sampling
    )

    c = RequestFetchCDSClient()

    try:
        with open(path / DATA_REQUESTS_FILENAME, "r") as fh:
            download_requests = yaml.load(fh, Loader=yaml.FullLoader)
    except FileNotFoundError:
        download_requests = {}

    for output_filename, query_kwargs in dl_queries:
        file_path = Path(path) / output_filename
        query_hash = utils.dict_to_hash(query_kwargs)

        download_id = str(file_path)
        should_make_request = True
        if not overwrite_existing and file_path.exists():
            if _data_valid(file_path=file_path, query_hash=query_hash):
                should_make_request = False
            else:
                warnings.warn(
                    "Invalid data found for ({})"
                    ", deleting and queuing for re-download"
                    "".format(output_filename)
                )
                file_path.unlink()

        if download_id in download_requests:
            if download_requests[download_id]["query_hash"] == query_hash:
                should_make_request = False
            else:
                del download_data[download_id]

        if should_make_request:
            request_id = c.queue_data_request(
                repository_name=REPOSITORY_NAME, query_kwargs=query_kwargs
            )
            assert request_id is not None
            download_requests[download_id] = dict(
                request_id=request_id, query_hash=query_hash
            )

    # save the download requests IDs so we have them later in case we quit,
    # downloads fail, etc
    if len(download_requests) > 0:
        Path(path).mkdir(exist_ok=True, parents=True)
        with open(path / DATA_REQUESTS_FILENAME, "w") as fh:
            fh.write(yaml.dump(download_requests))

    files_to_download = _get_files_to_download(path=path, c=c, debug=True)

    if len(files_to_download) > 0:
        print("Downloading files which are ready...")
        for file_path in files_to_download:
            request_details = download_requests[file_path]
            request_id = request_details["request_id"]
            query_hash = request_details["query_hash"]

            Path(file_path).parent.mkdir(exist_ok=True, parents=True)
            c.download_data_by_request(request_id=request_id, target=file_path)
            _fingerprint_downloaded_file(query_hash=query_hash, file_path=file_path)

            del download_requests[file_path]

    # save download requets again now we've downloaded the files that were
    # ready
    if len(download_requests) > 0:
        with open(path / DATA_REQUESTS_FILENAME, "w") as fh:
            fh.write(yaml.dump(download_requests))
    else:
        data_requests_file = Path(path / DATA_REQUESTS_FILENAME)
        if data_requests_file.exists():
            data_requests_file.unlink()
        print("All files downloaded!")


def all_data_is_downloaded(path):
    c = RequestFetchCDSClient()
    return len(_get_files_to_download(path=path, c=c)) == 0


def _get_files_to_download(path, c, debug=False):
    with open(path / DATA_REQUESTS_FILENAME, "r") as fh:
        download_requests = yaml.load(fh, Loader=yaml.FullLoader)

    files_to_download = []
    if len(download_requests) > 0:
        if debug:
            print("Status on current data requests:")
        for file_path, request_details in download_requests.items():
            request_id = request_details["request_id"]
            status = c.get_request_status(request_id=request_id)
            if debug:
                print(" {}:\n\t{} ({})".format(file_path, status, request_id))

            if status == "completed":
                files_to_download.append(file_path)
    return files_to_download


def _data_valid(file_path, query_hash):
    """
    Create hash on `query_kwargs` and ensure it matches that stored in file
    on disk
    """
    ds = netCDF4.Dataset(file_path)
    try:
        hash_unchanged = ds.getncattr("dict_checksum") == query_hash
    except AttributeError:
        # if hash checking fails, just download
        return False
    return hash_unchanged


def _build_query_times(model_run_type, t_start, t_end):
    """
    We currently download era5 data by date, build the query dates as datetime
    objects in the request time range
    Note: forecast needs to start a day earlier, due to it starting from 18Z
    """
    if model_run_type == "an":
        return [
            (t_start + datetime.timedelta(days=x))
            for x in range(0, (t_end - t_start).days + 1)
        ]
    elif model_run_type == "fc":
        return [
            (t_start + datetime.timedelta(days=x))
            for x in range(-1, (t_end - t_start).days + 1)
        ]
    else:
        raise NotImplementedError(model_run_type)


def _build_queries(t_start, t_end, bbox, latlon_sampling):
    """
    Return a list of cdsapi query arguments for any files that don't already
    exist
    """
    model_run_types = ["an", "fc"]  # analysis and forecast runs
    level_types = ["model", "single"]  # need model and surface data

    for model_run_type in model_run_types:
        query_times = _build_query_times(
            model_run_type=model_run_type, t_start=t_start, t_end=t_end
        )
        for level_type in level_types:
            for t in query_times:
                output_filename = FILENAME_FORMAT.format(
                    model_run_type=model_run_type,
                    level_type=level_type,
                    date=t.strftime(DATE_FORMAT),
                )

                query_kwargs = _build_query(
                    model_run_type=model_run_type,
                    level_type=level_type,
                    date=t,
                    bbox=bbox,
                    latlon_sampling=latlon_sampling,
                )

                yield (output_filename, query_kwargs)


def _build_query(model_run_type, level_type, date, bbox, latlon_sampling):
    if model_run_type == "an" and level_type == "single":
        return _build_single_level_an_query(
            date=date, bbox=bbox, latlon_sampling=latlon_sampling
        )
    elif model_run_type == "an" and level_type == "model":
        return _build_model_level_an_query(
            date=date, bbox=bbox, latlon_sampling=latlon_sampling
        )
    elif model_run_type == "fc" and level_type == "single":
        return _build_single_level_fc_query(
            date=date, bbox=bbox, latlon_sampling=latlon_sampling
        )
    elif model_run_type == "fc" and level_type == "model":
        return _build_model_level_fc_query(
            date=date, bbox=bbox, latlon_sampling=latlon_sampling
        )
    else:
        raise NotImplementedError(model_run_type, level_type)


def _get_query_status(request_id):
    return "TOFIX"


def _fingerprint_downloaded_file(query_hash, file_path):
    ds = netCDF4.Dataset(file_path, "a")
    ds.setncattr("dict_checksum", query_hash)
    ds.close()


def _build_single_level_an_query(date, bbox, latlon_sampling):
    # 2D ANALYSIS FROM ERA5
    # 31 Sea ice area fraction [(0 - 1)], ci
    # 32 Snow albedo [(0 - 1)], asn
    # 33 Snow density [kg m**-3], rsn
    # 34 Sea surface temperature [K], sst
    # 35 Ice temperature layer 1 [K], istl1
    # 39 Volumetric soil water layer 11 [m**3 m**-3], swvl1
    # 40 Volumetric soil water layer 21 [m**3 m**-3], swvl2
    # 41 Volumetric soil water layer 31 [m**3 m**-3], swvl3
    # 42 Volumetric soil water layer 41 [m**3 m**-3], swvl4
    # 129 Geopotential [m**2 s**-2], z
    # 134 Surface pressure [Pa], sp
    # 136 Total column water  [kg m**-2], tcw
    # 139 Soil temperature level 11 [K], stl1
    # 141 Snow depth [m of water equivalent], sd
    # 151 Mean sea level pressure [Pa ], msl
    # 159 Boundary layer height [m],  blh
    # 164 Total cloud cover [(0 - 1)],  tcc
    # 170 Soil temperature level 21 [K], stl2
    # 172 Land-sea mask [(0 - 1)], lsm
    # 183 Soil temperature level 31 [K], stl3
    # 186 Low cloud cover [(0 - 1)],  lcc
    # 187 Medium cloud cover  [(0 - 1)],  mcc
    # 188 High cloud cover  [(0 - 1)],  hcc
    # 236 Soil temperature level 41 [K], stl4
    # 238 Temperature of snow layer [K], tsn
    # 243 Forecast albedo [(0 - 1)], fal
    # 244 Forecast surface roughness [m], fsr
    # 245 Forecast logarithm of surface roughness for heat [~], flsr

    return {
        "class": "ea",
        "date": date.strftime(DATE_FORMAT),
        "expver": "1",
        "levtype": "sfc",
        "param": (
            "31.128/32.128/33.128/34.128/35.128/39.128/40.128/41.128"
            "/42.128/129.128/136.128/134.128/139.128/141.128/151.128"
            "/159.128/164.128/170.128/172.128/183.128/186.128/187.128"
            "/188.128/236.128/238.128/243.128/244.128/245.128"
        ),
        "stream": "oper",
        "time": (
            "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00"
            "/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00"
            "/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00"
            "/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00"
        ),
        "area": [bbox.lat_max, bbox.lon_min, bbox.lat_min, bbox.lon_max,],
        "grid": "{}/{}".format(latlon_sampling.lat, latlon_sampling.lon),
        "type": "an",
        "format": "netcdf",
    }


def _build_model_level_an_query(date, bbox, latlon_sampling):
    # 3D ANALYSIS FROM ERA5
    # 75  Specific rain water content [kg kg**-1],  crwc
    # 76  Specific snow water content [kg kg**-1],  cswc
    # 77  Eta-coordinate vertical velocity [s**-1], etadot
    # 129 Geopotential* [m**2 s**-2], z
    # 130 Temperature [K],  t
    # 131 U component of wind [m s**-1],  u
    # 132 V component of wind [m s**-1],  v
    # 133 Specific humidity [kg kg**-1],  q
    # 135 Vertical velocity [Pa s**-1], w
    # 152 Logarithm of surface pressure*  [~],  lnsp
    # 203 Ozone mass mixing ratio [kg kg**-1],  o3
    # 246 Specific cloud liquid water content [kg kg**-1],  clwc
    # 247 Specific cloud ice water content  [kg kg**-1],  ciwc
    # 248 Fraction of cloud cover [(0 - 1)],  cc
    return {
        "class": "ea",
        "date": date.strftime(DATE_FORMAT),
        "expver": "1",
        "levelist": (
            "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/"
            "22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/"
            "40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/"
            "58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/"
            "76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/"
            "94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/"
            "109/110/111/112/113/114/115/116/117/118/119/120/121/"
            "122/123/124/125/126/127/128/129/130/131/132/133/134/"
            "135/136/137"
        ),
        "levtype": "ml",
        "param": "75/76/77/129/130/131/132/133/135/152/203/246/247/248",
        "stream": "oper",
        "time": (
            "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00"
            "/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00"
            "/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00"
            "/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00"
        ),
        "area": [bbox.lat_max, bbox.lon_min, bbox.lat_min, bbox.lon_max],
        "grid": "{}/{}".format(latlon_sampling.lat, latlon_sampling.lon),
        "type": "an",
        "format": "netcdf",
    }


def _build_single_level_fc_query(date, bbox, latlon_sampling):
    # 2D FORECAST DATA FROM ERA5
    # 228001 Convective inhibition [J kg**-1], cin
    # 228003 Friction velocity [m s**-1], zust
    # 228023 Cloud base height [m], cbh
    # 235033 Mean surface sensible heat flux [W m**-2], msshf
    # 235034 Mean surface latent heat flux [W m**-2], mslhf
    # 235035 Mean surface downward short-wave radiation flux [W m**-2],
    #        msdwswrf
    # 235036 Mean surface downward long-wave radiation flux [W m**-2], msdwlwrf
    # 235037 Mean surface net short-wave radiation flux [W m**-2], msnswrf
    # 235038 Mean surface net long-wave radiation flux [W m**-2], msnlwrf
    # 235039 Mean top net short-wave radiation flux [W m**-2], mtnswrf
    # 235040 Mean top net long-wave radiation flux [W m**-2], mtnlwrf
    # 235043 Mean evaporation rate [kg m**-2 s**-1 ], mer
    # 235049 Mean top net short-wave radiation flux, clear sky [W m**-2],
    #        mtnswrfcs
    # 235050 Mean top net long-wave radiation flux, clear sky [W m**-2],
    #        mtnlwrfcs
    # 235051 Mean surface net short-wave radiation flux, clear sky [W m**-2],
    #        msnswrfcs
    # 235052 Mean surface net long-wave radiation flux, clear sky [W m**-2],
    #        msnlwrfcs
    # 235053 Mean top downward short-wave radiation flux [W m**-2], mtdwswrf
    # 235058 Mean surface direct short-wave radiation flux [W m**-2], msdrswrf
    # 235059 Mean surface direct short-wave radiation flux, clear sky [W m**-2]
    #        msdrswrfcs
    # 235068 Mean surface downward short-wave radiation flux, clear sky
    #        [W m**-2],msdwswrfcs
    # 235069 Mean surface downward long-wave radiation flux, clear sky
    #        [W m**-2], msdwlwrfcs
    # 235070 Mean potential evaporation rate [kg m**-2 s**-1 ], mper
    return {
        "class": "ea",
        "date": date.strftime(DATE_FORMAT),
        "expver": "1",
        "levtype": "sfc",
        "param": (
            "228001/228003/228023/235033/235034/235035/235036/235037/"
            "235038/235039/235040/235043/235049/235050/235051/235052/"
            "235053/235058/235059/235068/235069/235070"
        ),
        "stream": "oper",
        "time": "06:00:00/18:00:00",
        "area": [bbox.lat_max, bbox.lon_min, bbox.lat_min, bbox.lon_max],
        "grid": "{}/{}".format(latlon_sampling.lat, latlon_sampling.lon),
        "type": "fc",
        "step": "0/1/2/3/4/5/6/7/8/9/10/11",
        "format": "netcdf",
    }


def _build_model_level_fc_query(date, bbox, latlon_sampling):
    # 3D FORECAST DATA FROM ERA5
    # 235001  Mean temperature tendency due to short-wave radiation [K s**-1],
    #         mttswr
    # 235002  Mean temperature tendency due to long-wave radiation  [K s**-1],
    #         mttlwr
    # 235003  Mean temperature tendency due to short-wave radiation
    #         [clear sky  K s**-1], mttswrcs
    # 235004  Mean temperature tendency due to long-wave radiation
    #         [clear sky K s**-1], mttlwrcs
    return {
        "class": "ea",
        "date": date.strftime(DATE_FORMAT),
        "expver": "1",
        "levelist": (
            "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/"
            "22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/"
            "40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/"
            "58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/"
            "76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/"
            "94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/"
            "109/110/111/112/113/114/115/116/117/118/119/120/121/"
            "122/123/124/125/126/127/128/129/130/131/132/133/134/"
            "135/136/137"
        ),
        "levtype": "ml",
        "param": "235001/235002/235003/235004",
        "stream": "oper",
        "time": "06:00:00/18:00:00",
        "type": "fc",
        "time": "06:00:00/18:00:00",
        "area": [bbox.lat_max, bbox.lon_min, bbox.lat_min, bbox.lon_max,],
        "grid": "{}/{}".format(latlon_sampling.lat, latlon_sampling.lon),
        "step": "0/1/2/3/4/5/6/7/8/9/10/11",
        "format": "netcdf",
    }


def _era_5_normalise_longitude(ds):
    """Normalise longitudes to be between 0 and 360 degrees
    This is needed because these are stored differently in the surface
    and model level data. Rounding up to 4 decimals seems to work for now,
    with more decimals misalignment has happenend. Would be good to sort
    out why this is the case.
    """

    def longitude_set_meridian(longitude):
        """Sets longitude to be between -180 and 180 degrees"""
        return (longitude + 180.0) % 360.0 - 180.0

    ds.coords["longitude"] = (
        "longitude",
        np.round(longitude_set_meridian(ds.coords["longitude"]), decimals=4),
        ds.coords["longitude"].attrs,
    )
    return ds


def load_data(data_path):
    datasets = {}

    model_run_types = ["an", "fc"]  # analysis and forecast runs
    level_types = ["model", "single"]  # need model and surface data

    for model_run_type in model_run_types:
        datasets_run = []
        for level_type in level_types:
            filename_format = FILENAME_FORMAT.format(
                model_run_type=model_run_type, level_type=level_type, date="*"
            )

            files = data_path.glob(filename_format)

            ds_ = xr.open_mfdataset(files, combine="by_coords")
            # z needs to be dropped to prevent duplicity, lnsp is simply
            # redundant
            if model_run_type == "an" and level_type == "model":
                ds_ = ds_.drop_vars(["lnsp"])
            ds_ = _era_5_normalise_longitude(ds=ds_)
            datasets_run.append(ds_)
        ds_run = xr.merge(datasets_run, compat="override")
        datasets[model_run_type] = ds_run

    ds = xr.merge(datasets, compat="override")
    ds = ds.rename(dict(latitude="lat", longitude="lon"))

    return ds
