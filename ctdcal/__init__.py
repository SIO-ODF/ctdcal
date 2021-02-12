"""
Library and scripts for processing CTD profiles with bottle
data, and similar data types such as ARGO and moored CTD
sensors.
"""

from . import *
from ._version import get_versions
import os.path
import types

__version__ = get_versions()["version"]
del get_versions


def get_ctdcal_config():
    """
    Find and load config file.

    Directory search mimicked from cchdo/hdo-uow
    File read mimicked from pallets/flask
    """
    last_dir = os.getcwd()
    filepath = None
    config_file = "config.py"

    if os.path.exists(os.path.join(last_dir, config_file)):
        filepath = os.path.join(last_dir, config_file)
    else:
        while last_dir != os.path.dirname(last_dir):
            last_dir = os.path.dirname(last_dir)
            if os.path.exists(os.path.join(last_dir, config_file)):
                filepath = os.path.join(last_dir, config_file)
                break
    if not filepath:
        raise FileNotFoundError(f"Failed to find config.py file in {os.getcwd()}")

    # compile config.py and save variables to dict
    config = {}
    try:
        with open(filepath, mode="rb") as f:
            exec(compile(f.read(), filepath, "exec"), config)
    except OSError:
        print("Failed to load config file")

    for k in list(config.keys()):
        if k.startswith("__"):
            del config[k]

    return type("config", (object,), config)  # make into class for easier calls
