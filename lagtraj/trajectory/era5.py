import yaml
import datetime
import xarray as xr
import os
import numpy as np
import pandas as pd
from lagtraj.utils.parsers import domain_filename_parse,trajectory_filename_parse

def main():
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('input_file')
    argparser.add_argument ("-f", "--file", dest='directories_file', default='directories.yaml', type=str)
    argparser.add_argument ("-o", "--overwrite", dest='l_overwrite',action='store_true')
    args = argparser.parse_args()
    get_from_yaml(args.input_file,args.directories_file,args.l_overwrite)

def get_from_yaml(input_file,directories_file,l_overwrite):
    with open(directories_file) as this_directories_file:
        directories_dict=yaml.load(this_directories_file, Loader=yaml.FullLoader)  
    with open(input_file) as this_trajectory_file:
        trajectory_dict=yaml.load(this_trajectory_file, Loader=yaml.FullLoader)
    
    trajectory_type=(trajectory_dict['trajectory_type']).lower()
    if(trajectory_type=='eulerian'):
        create_eulerian_trajectory(directories_dict,trajectory_dict)
    elif(trajectory_type=='single_level'):
        create_single_level_trajectory(directories_dict,trajectory_dict)
    elif(trajectory_type=='weighted'):
        create_weighted_trajectory(directories_dict,trajectory_dict)
    else:
        raise ValueError('Trajectory_type not found')

def create_eulerian_trajectory(directories_dict,trajectory_dict):
    time_list=pd.date_range(np.datetime64(trajectory_dict['datetime_end'])-np.timedelta64(24, 'h'), periods=25, freq='h').index.values.astype(float)
    nr_hours=len(time_list)
    lats=np.full((nr_hours),trajectory_dict['lat_end'])
    lons=np.full((nr_hours),trajectory_dict['lon_end']%360)
    attrs = {'units': 'hours since 2000-01-01'}
    data=xr.Dataset({'time': ('time', time_list, attrs)})
    data['lat'] = (('time'), lats)                 
    data['lon'] = (('time'), lons)  
    var_attrs = {
        'lon': {'long_name': 'longitude',
                'units': 'degrees_north'},
        'lat': {'long_name': 'latitude',
                'units': 'degrees_east'},
         }
    for this_var,this_attr in var_attrs.items():
        data[this_var] = data[this_var].assign_attrs(**this_attr)
    data.to_netcdf(trajectory_filename_parse(directories_dict,trajectory_dict))
      
def create_single_level_trajectory(directories_dict,trajectory_dict):
    pass

def create_weighted_trajectory(directories_dict,trajectory_dict):
    pass

def trajectory_to_xarray():
    pass
    
if __name__ == '__main__':
    main()
