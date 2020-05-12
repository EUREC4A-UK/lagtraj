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
    args = argparser.parse_args()
    get_from_yaml(args.input_file,args.directories_file,args.l_overwrite)

def get_from_yaml(input_file,directories_file,l_overwrite):
    with open(directories_file) as this_directories_file:
        directories_dict=yaml.load(this_directories_file, Loader=yaml.FullLoader)  
    with open(input_file) as this_forcing_file:
        forcing_dict=yaml.load(this_forcing_file, Loader=yaml.FullLoader)
    
    
    
if __name__ == '__main__':
    main()
