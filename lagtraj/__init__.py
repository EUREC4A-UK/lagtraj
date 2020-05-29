from pathlib import Path
import os


# by default we store data relative to where lagtraj is invoked from
DEFAULT_DATA_PATH = Path(os.getcwd())/"data"
