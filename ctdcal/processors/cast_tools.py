"""
Tools for cleaning up converted and processed cast data.
"""
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import signal as sig


class Cast(object):
    """
    Cast data container with methods for separating upcast, filtering
    and trimming soak period.
    """
    def __init__(self, cast_id, datadir):
        self.cast_id = cast_id  # cast identifier
        self.datadir = datadir  # data directory
        self.p_col = None  # name of pressure column
        self.proc = None  # processed cast data
        self.downcast = None
        self.filtered = None
        self.trimmed = None
        self.load_cast()

    def load_cast(self):
        """
        Read the processed data into a dataframe.
        """
        f = Path(self.datadir, 'processed/%s.pkl' % self.cast_id)
        self.proc = pd.read_pickle(f)

    def parse_downcast(self, data):
        """
        Trim away the upcast starting from the max value of the pressure
        column.

        Parameters
        ----------
        data : DataFrame
            cast data
        """
        if self.p_col is None:
            raise AttributeError("Pressure column (p_col) attribute is not set")
        downcast = data.loc[:self.proc[self.p_col].idxmax()].copy()
        self.downcast = downcast

    def filter(self, data, win_size, win_type="triangle", cols=None):
        """
        Filter processed CTD data using one of three window types (boxcar,
        hann, triangle).

        Parameters
        ----------
        data : DataFrame, array-like
            cast data
        win_size : int
            window size in number of samples
        win_type : str
            window type
        cols : list
            column names to filter
        """
        data = data.copy()
        if win_type == "boxcar":
            win = sig.windows.boxcar(win_size)
        elif win_type == "hann":
            win = sig.windows.hann(win_size)
        elif win_type == "triangle":
            win = sig.windows.triang(win_size)
        else:
            raise AttributeError('Error filtering cast %s! No filter of type %s found.' % (self.cast_id, win_type))
        for col in cols:
            data[col] = sig.convolve(data[col], win, mode="same") / np.sum(win)
        self.filtered = data

    def trim_soak(self, data, win_size, max_soak=20):
        """
        Detects the end of the soak period by identifying local peaks in the
        cast data, and choosing the last one above the threshold.

        data : DataFrame
            cast data
        win_size : int
            window size in number of samples
        max_soak : int
            soak threshold in cast pressure units
        """
        if self.p_col is None:
            raise AttributeError("Pressure column (p_col) attribute is not set")
        if self.downcast is None:
            raise AttributeError("Downcast attribute is not set.")
        if win_size < 0:
            win_size = 0  # n must be positive or zero

        # Adapted from: https://stackoverflow.com/questions/48023982/pandas-finding-local-max-and-min
        # Find local peaks
        local_min_vals = data.loc[data[self.p_col] == data[self.p_col].rolling(win_size, center=True).min(), [self.p_col]]

        # apply max_soak threshold to filter false positives...
        i = local_min_vals.index[-1]
        while local_min_vals.loc[i, self.p_col] > max_soak:
            try:
                local_min_vals.drop(i, inplace=True)
            except KeyError:
                print('Whoa, trouble finding the soak on cast %s! Nothing was trimmed!' % self.cast_id)
                return
            i = local_min_vals.index[-1]

        soak_end = local_min_vals.index[-1]
        self.trimmed = self.downcast.loc[soak_end:]
