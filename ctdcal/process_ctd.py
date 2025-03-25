"""
A module for handling continuous CTD data processing, including file write-outs.
"""

import logging
import warnings
from datetime import datetime, timezone
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from . import get_ctdcal_config, oxy_fitting
from .processors.functions_oxy import calculate_dV_dt

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

warnings.filterwarnings("ignore", "Mean of empty slice.")


