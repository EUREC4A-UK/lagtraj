import glob
import importlib
import os


def _package_contents():
    for path in glob.glob(os.path.join(os.path.dirname(__file__), "*.py")):
        path = os.path.basename(path)
        if not path.startswith("_"):
            module_name = path.replace(".py", "")
            yield module_name, importlib.import_module(f"{__package__}.{module_name}")


available = dict(_package_contents())
