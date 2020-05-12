import yaml
import datetime
import xarray as xr
import os
import numpy as np
import pandas as pd
from lagtraj.utils.parsers import domain_filename_parse,trajectory_filename_parse,forcings_filename_parse

def main():
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('input_file')
    argparser.add_argument ("-f", "--file", dest='directories_file', default='directories.yaml', type=str)
    args = argparser.parse_args()
    get_from_yaml(args.input_file,args.directories_file)

def get_from_yaml(input_file,directories_file):
    with open(directories_file) as this_directories_file:
        directories_dict=yaml.load(this_directories_file, Loader=yaml.FullLoader)  
    with open(input_file) as this_forcings_file:
        forcings_dict=yaml.load(this_forcings_file, Loader=yaml.FullLoader)  
    # GET TRAJECTORY FILE
    
    # FOR EACH TIME STEP
    # reinterpolate data to height or pressure with effective height (tbd)
    # create 'mask' for forcings on domain
    # calculate profiles and forcings
    
def append_timestep():
    pass

def export_to_hightune():
    pass
        
if __name__ == '__main__':
    main()
