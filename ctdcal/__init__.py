"""
Library and scripts for processing CTD profiles with bottle
data, and similar data types such as ARGO and moored CTD
sensors.
"""

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from . import *
