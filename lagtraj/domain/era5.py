import cdsapi
from os import path
import yaml
#from utils import 
cds_client = cdsapi.Client()

def main():
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('input_file')
    argparser.add_argument('start_date')
    argparser.add_argument('end_date')
    argparser.add_argument ("-f", "--file", dest='directories_file', default='directories.yaml', type=str);
    args = argparser.parse_args()
    get_from_yaml(args.input_file,args.start_date,args.end_date,args.directories_file)

def get_from_yaml(input_file):
    with open(input_file) as this_case_file:
        domain_dict=yaml.load(this_case_file, Loader=yaml.FullLoader)
    print(domain_dict)
    #for date in dates:
        #~ get_by_day(date,north,west,south,east,era5dir)
        
def get_by_date(date,north,west,south,east,era5dir,overwrite=False):
    date_exists=check_if_date_exists(date,north,west,south,east,era5dir)
    if((overwrite==True) or (not day_exists)):
        get_single_level_data(datestr,north,west,south,east)
        get_model_level_data(datestr,north,west,south,east)
   
def get_model_level_data(date,north,west,south,east):
    c.retrieve('reanalysis-era5-complete', {
    'class': 'ea',
    'date': datestr,
    'expver': '1',
    'levelist': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
    'levtype': 'ml',
    'param': '75/76/77/129/130/131/132/133/135/152/203/246/247/248',
    'stream': 'oper',
    'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
    'type': 'an',
    'step': '0',
    'area': [north, west, south, east],
    'format':'netcdf',
    }, era5dir+'/ml_'+date+'_'+str(north)+'_'+str(west)+'_'+str(south)+'_'+str(east)+'.nc')
    
def get_single_level_data():
    pass

def check_if_day_exists(date,north,west,south,east,era5dir):
    pass
    #domain_filename_parse
    #path.isfile()

def domain_filename_parse(date,north,west,south,east,era5dir):
    filename='domain_'+date
    return filename

if __name__ == '__main__':
    main()
