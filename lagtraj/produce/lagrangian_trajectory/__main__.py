import dateutil.parser
import argparse


from . import main


argparser = argparse.ArgumentParser()
argparser.add_argument("lat0", help="starting latitude")
argparser.add_argument("lon0", help="starting longitude")
argparser.add_argument("t0", help="starting time", type=dateutil.parser.parse)
argparser.add_argument("t_max", help="end time", type=dateutil.parser.parse)
argparser.add_argument("dt", help="time increment")
argparser.add_argument("--out-filename", help="output filename (netCDF)")

args = argparser.parse_args()

main(**dict(args))
