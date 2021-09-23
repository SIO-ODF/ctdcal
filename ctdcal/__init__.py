"""
Library and scripts for processing CTD profiles with bottle
data, and similar data types such as ARGO and moored CTD
sensors.
"""

import logging
from importlib import resources
from importlib.metadata import version, PackageNotFoundError

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
    with resources.path("ctdcal", "config.py") as filepath:
        try:
            with open(filepath, mode="rb") as f:
                exec(compile(f.read(), filepath, "exec"), config)
        except OSError:
            log.error(f"Failed to load config file {filepath}")

    for k in list(config.keys()):
        if k.startswith("__"):
            del config[k]

    return type("config", (object,), config)  # make into class for easier calls
