"""
Module for processing oxygen from CTD and bottle samples.
"""

import csv
import logging
import xml.etree.cElementTree as ET
from collections import OrderedDict
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import scipy

from ctdcal.fitting.fit_common import NodeNotFoundError, get_node

from . import flagging as flagging
from . import get_ctdcal_config
from . import process_ctd as process_ctd
from ctdcal.plotting.plot_fit import _intermediate_residual_plot






