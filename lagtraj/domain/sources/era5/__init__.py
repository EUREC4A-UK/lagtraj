FILENAME_FORMAT = "{model_run_type}_{level_type}_{date}.nc"

VERSION_FILENAME = "VERSION"


from .download import download_data, all_data_is_downloaded  # noqa
from .load import load_data  # noqa
