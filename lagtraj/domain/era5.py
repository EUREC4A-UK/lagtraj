import cdsapi
import multiprocessing
import yaml
from lagtraj.utils.parsers import domain_filename_parse
import datetime
import netCDF4
import os
import zlib

# Routines for downloading data from ECMWF archives
# TODO
# - Think about a way to put multiple days in one file, which may speed up downloading?
#   (on request: check for each day individually if it is there already, add the missing ones to a list, download the request, split the files by date again) 
# - Make both NetCDF and Grib downloads are an option? 

cds_client = cdsapi.Client()

def main():
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('input_file')
    argparser.add_argument('start_date')
    argparser.add_argument('end_date')
    argparser.add_argument ("-f", "--file", dest='directories_file', default='directories.yaml', type=str)
    argparser.add_argument ("-o", "--overwrite", dest='l_overwrite',action='store_true', default=False)
    args = argparser.parse_args()
    get_from_yaml(args.input_file,args.start_date,args.end_date,args.directories_file,args.l_overwrite)

def get_from_yaml(input_file,start_date,end_date,directories_file,l_overwrite):
    with open(directories_file) as this_directories_file:
        directories_dict=yaml.load(this_directories_file, Loader=yaml.FullLoader)  
    with open(input_file) as this_domain_file:
        domain_dict=yaml.load(this_domain_file, Loader=yaml.FullLoader)
            
    start = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    end = datetime.datetime.strptime(end_date, "%Y-%m-%d")
    # Note: forecast needs to start a day earlier, due to it starting from 18Z
    an_dates_as_list = [(start + datetime.timedelta(days=x)).strftime("%Y-%m-%d") for x in range(0, (end-start).days+1)]
    fc_dates_as_list = [(start + datetime.timedelta(days=x)).strftime("%Y-%m-%d") for x in range(-1, (end-start).days+1)]

    for this_date in an_dates_as_list:
        get_single_level_an(this_date,directories_dict,domain_dict,l_overwrite)
        get_model_level_an(this_date,directories_dict,domain_dict,l_overwrite)

    for this_date in fc_dates_as_list:
        get_single_level_fc(this_date,directories_dict,domain_dict,l_overwrite)
        get_model_level_fc(this_date,directories_dict,domain_dict,l_overwrite)


def create_hash(this_dict):
    this_hash=0
    for item in sorted(this_dict.items()):
        curr_hash = 1
        for sub_item in item:
            curr_hash = zlib.adler32(bytes(repr(sub_item),'utf-8'), curr_hash)
        this_hash=this_hash ^ curr_hash
    return str(this_hash)
    
def retrieve_file(this_repository,this_dict,this_filename):
    this_hash=create_hash(this_dict)
    c = cdsapi.Client()
    c.retrieve(this_repository,this_dict,this_filename)
    ds=netCDF4.Dataset(this_filename,'a')
    ds.setncattr('dict_checksum', this_hash)
    ds.close()
    
def get_single_level_an(date,directories_dict,domain_dict,l_overwrite):        
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
    prefix_str='single_an'
    this_filename=domain_filename_parse(date,prefix_str,directories_dict,domain_dict)
    this_repository='reanalysis-era5-complete'
    this_dict={
        'class': 'ea',
        'date': date,
        'expver': '1',
        'levtype': 'sfc',
        'param': '31.128/32.128/33.128/34.128/35.128/39.128/40.128/41.128/42.128/129.128/136.128/134.128/139.128/141.128/151.128/159.128/164.128/170.128/172.128/183.128/186.128/187.128/188.128/236.128/238.128/243.128/244.128/245.128',
        'stream': 'oper',
        'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
        'area': [domain_dict['lat_max'], domain_dict['lon_min'], domain_dict['lat_min'], domain_dict['lon_max']],
        'grid': str(domain_dict['lat_samp'])+'/'+str(domain_dict['lon_samp']),
        'type': 'an',
        'format':'netcdf'}
    file_exists=check_if_file_exists(date,prefix_str,directories_dict,domain_dict,this_dict)
    if(l_overwrite or not(file_exists)):       
        p = multiprocessing.Process(target=retrieve_file, args=(this_repository,this_dict,this_filename,))
        p.start()
    
def get_model_level_an(date,directories_dict,domain_dict,l_overwrite):
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
    prefix_str='model_an'
    this_filename=domain_filename_parse(date,prefix_str,directories_dict,domain_dict)
    this_repository='reanalysis-era5-complete'
    this_dict={
        'class': 'ea',
        'date': date,
        'expver': '1',
        'levelist': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
        'levtype': 'ml',
        'param': '75/76/77/129/130/131/132/133/135/152/203/246/247/248',
        'stream': 'oper',
        'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
        'area': [domain_dict['lat_max'], domain_dict['lon_min'], domain_dict['lat_min'], domain_dict['lon_max']],
        'grid': str(domain_dict['lat_samp'])+'/'+str(domain_dict['lon_samp']),
        'type': 'an',
        'format':'netcdf'}
    file_exists=check_if_file_exists(date,prefix_str,directories_dict,domain_dict,this_dict)
    if(l_overwrite or not(file_exists)):       
        p = multiprocessing.Process(target=retrieve_file, args=(this_repository,this_dict,this_filename,))
        p.start()

def get_single_level_fc(date,directories_dict,domain_dict,l_overwrite):
    # 2D FORECAST DATA FROM ERA5
    # 228001 Convective inhibition [J kg**-1], cin
    # 228003 Friction velocity [m s**-1], zust
    # 228023 Cloud base height [m], cbh
    # 235033 Mean surface sensible heat flux [W m**-2], msshf
    # 235034 Mean surface latent heat flux [W m**-2], mslhf
    # 235035 Mean surface downward short-wave radiation flux [W m**-2], msdwswrf
    # 235036 Mean surface downward long-wave radiation flux [W m**-2], msdwlwrf
    # 235037 Mean surface net short-wave radiation flux [W m**-2], msnswrf
    # 235038 Mean surface net long-wave radiation flux [W m**-2], msnlwrf
    # 235039 Mean top net short-wave radiation flux [W m**-2], mtnswrf
    # 235040 Mean top net long-wave radiation flux [W m**-2], mtnlwrf
    # 235043 Mean evaporation rate [kg m**-2 s**-1 ], mer
    # 235049 Mean top net short-wave radiation flux, clear sky [W m**-2], mtnswrfcs
    # 235050 Mean top net long-wave radiation flux, clear sky [W m**-2], mtnlwrfcs
    # 235051 Mean surface net short-wave radiation flux, clear sky [W m**-2], msnswrfcs
    # 235052 Mean surface net long-wave radiation flux, clear sky [W m**-2], msnlwrfcs
    # 235053 Mean top downward short-wave radiation flux [W m**-2], mtdwswrf
    # 235058 Mean surface direct short-wave radiation flux [W m**-2], msdrswrf
    # 235059 Mean surface direct short-wave radiation flux, clear sky [W m**-2], msdrswrfcs
    # 235068 Mean surface downward short-wave radiation flux, clear sky [W m**-2], msdwswrfcs
    # 235069 Mean surface downward long-wave radiation flux, clear sky [W m**-2], msdwlwrfcs
    # 235070 Mean potential evaporation rate [kg m**-2 s**-1 ], mper
    prefix_str='single_fc'
    this_filename=domain_filename_parse(date,prefix_str,directories_dict,domain_dict)
    this_repository='reanalysis-era5-complete'
    this_dict={
        'class': 'ea',
        'date': date,
        'expver': '1',
        'levtype': 'sfc',
        'param': '228001/228003/228023/235033/235034/235035/235036/235037/235038/235039/235040/235043/235049/235050/235051/235052/235053/235058/235059/235068/235069/235070',
        'stream': 'oper',
        'time': '06:00:00/18:00:00',
        'area': [domain_dict['lat_max'], domain_dict['lon_min'], domain_dict['lat_min'], domain_dict['lon_max']],
        'grid': str(domain_dict['lat_samp'])+'/'+str(domain_dict['lon_samp']),
        'type': 'fc',
        'step':'0/1/2/3/4/5/6/7/8/9/10/11',
        'format':'netcdf'}
    file_exists=check_if_file_exists(date,prefix_str,directories_dict,domain_dict,this_dict)
    if(l_overwrite or not(file_exists)): 
        p = multiprocessing.Process(target=retrieve_file, args=(this_repository,this_dict,this_filename,))
        p.start()
            
def get_model_level_fc(date,directories_dict,domain_dict,l_overwrite):       
    # 3D FORECAST DATA FROM ERA5
    # 235001  Mean temperature tendency due to short-wave radiation [K s**-1],  mttswr
    # 235002  Mean temperature tendency due to long-wave radiation  [K s**-1],  mttlwr
    # 235003  Mean temperature tendency due to short-wave radiation [clear sky  K s**-1], mttswrcs
    # 235004  Mean temperature tendency due to long-wave radiation [clear sky K s**-1], mttlwrcs
    prefix_str='model_fc'
    this_filename=domain_filename_parse(date,prefix_str,directories_dict,domain_dict)
    this_repository='reanalysis-era5-complete'
    this_dict={
        'class': 'ea',
        'date': date,
        'expver': '1',
        'levelist': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
        'levtype': 'ml',
        'param': '235001/235002/235003/235004',
        'stream': 'oper',
        'time': '06:00:00/18:00:00',
        'type': 'fc',
        'time': '06:00:00/18:00:00',
        'area': [domain_dict['lat_max'], domain_dict['lon_min'], domain_dict['lat_min'], domain_dict['lon_max']],
        'grid': str(domain_dict['lat_samp'])+'/'+str(domain_dict['lon_samp']),
        'step':'0/1/2/3/4/5/6/7/8/9/10/11',
        'format':'netcdf'}
    file_exists=check_if_file_exists(date,prefix_str,directories_dict,domain_dict,this_dict)
    if(l_overwrite or not(file_exists)):       
        p = multiprocessing.Process(target=retrieve_file, args=(this_repository,this_dict,this_filename,))
        p.start()
        
def check_if_file_exists(date,prefix_str,directories_dict,domain_dict,this_dict):
    this_filename=domain_filename_parse(date,prefix_str,directories_dict,domain_dict)
    # check if the file name exists
    if(not(os.path.isfile(this_filename))):
        return False
    else:
        # check if the meta data agrees
        this_hash=create_hash(this_dict)
        ds=netCDF4.Dataset(this_filename)
        hash_unchanged=(ds.getncattr('dict_checksum')==this_hash)
        if(hash_unchanged):
            return True
        else:
            return False

if __name__ == '__main__':
    main()
