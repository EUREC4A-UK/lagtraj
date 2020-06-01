from pathlib import Path
import dateutil.parser


INPUT_REQUIRED_FIELDS = {
    'trajectory_type': str,
    'domain': str,
    'lat_origin': float,
    'lon_origin': float,
    'datetime_origin': dateutil.parser.parse,
    'forward_duration_hours|backward_duration_hours': int,
}


def build_data_path(data_path, trajectory_name):
    return Path(data_path)/"{}.nc".format(trajectory_name)
