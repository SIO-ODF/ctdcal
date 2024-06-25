"""
Library and scripts for processing CTD profiles with bottle
data, and similar data types such as ARGO and moored CTD
sensors.
"""

import logging
import pathlib
from importlib import resources
from importlib.metadata import PackageNotFoundError, version

log = logging.getLogger(__name__)

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    pass


def get_ctdcal_config():
    """
    Find and load config file from ctdcal module folder (ctdcal/ctdcal/).

    File read mimicked from pallets/flask
    """
    # compile config.py and save variables to dict
    config = {}
    resource_path = pathlib.Path(resources.files("ctdcal"))
    config_file_path = resource_path / "config.py"

    try:
        #   Read the config file as bytes, compile the bytes, and create the 'config' dictionary
        exec(compile(config_file_path.read_bytes(), str(config_file_path), "exec"), config)
    except OSError:
        log.error(f"Failed to load config file {config_file_path}")

    for k in list(config.keys()):
        if k.startswith("__"):
            del config[k]

    return type("config", (object,), config)  # make into class for easier calls
