#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:package: ctdcal.tools.cast_tools
:file: ctdcal/tools/cast_tools.py
:author: Allen Smith
:brief: Tools for cleaning up converted and processed cast data.
"""
from pathlib import Path

import pandas as pd
import numpy as np
from munch import munchify
from scipy import signal as sig

from ctdcal import get_ctdcal_config as cfg


# configs
dirs = munchify(cfg().dirs)
deltas = munchify(cfg().despike_deltas)

# constants
FLT_WIN_SIZE = cfg().win_size
MAX_SOAK = cfg().max_soak


class Cast(object):
    """
    Cast processing and cleaning class with methods for trimming, filtering
    and despiking.
    """
    def __init__(self, cast_num):
        self.cast_num = cast_num
        self.raw = None
        self.downcast = None
        self.filtered = None
        self.load_cast()

    def load_cast(self):
        """
        Read the converted data into a dataframe.
        """
        f = Path(dirs.converted, '%s.pkl' % self.cast_num)
        self.raw = pd.read_pickle(f)

    def parse_downcast(self, data):
        """
        Trim away the upcast starting from the max pressure value and add a datetime column.

        :param data: dataframe of cast
        """
        downcast = data.loc[:self.raw['CTDPRS'].idxmax()].copy()
        downcast['time'] = pd.to_datetime(downcast['scan_datetime'].dropna().astype(np.int64), unit='s')
        self.downcast = downcast

    def filter(self, data, window="triangle", win_size=FLT_WIN_SIZE, cols=None):
        """
        Filter raw CTD data using one of three window types (boxcar, hann, triangle).

        :param data: DataFrame, cast data
        :param window: str, filter window type
        :param win_size: int, window size in seconds
        :param cols: list, data columns to process
        """
        data = data.copy()
        win_size = win_size * 24
        for col in cols:
            if window == "boxcar":
                win = sig.windows.boxcar(win_size)
            elif window == "hann":
                win = sig.windows.hann(win_size)
            elif window == "triangle":
                win = sig.windows.triang(win_size)
            else:
                print('Error filtering cast %s! Filter "%s" not found' % (self.cast_num, window))
                break
            data[col] = sig.convolve(data[col], win, mode="same") / np.sum(win)
        self.filtered = data

    def trim_soak(self, data, n=20):
        """
        Detects the end of the soak period by identifying local peaks in the
        cast data, and choosing the last one above the threshold.

        :param data: DataFrame, cast data
        :param n: int, window size in seconds
        """
        n = n * 24
        if n < 0:
            n = 0  # n must be positive or zero

        # Adapted from: https://stackoverflow.com/questions/48023982/pandas-finding-local-max-and-min
        # Find local peaks
        local_min_vals = data.loc[data['CTDPRS'] == data['CTDPRS'].rolling(n, center=True).min(), ['CTDPRS']]

        # apply MAX_SOAK threshold to filter false positives...
        i = local_min_vals.index[-1]
        while local_min_vals.loc[i, 'CTDPRS'] > MAX_SOAK:
            try:
                local_min_vals.drop(i, inplace=True)
            except KeyError:
                print('Whoa, trouble finding the soak on cast %s! Nothing was trimmed!' % self.cast_num)
                return
            i = local_min_vals.index[-1]

        soak_end = local_min_vals.index[-1]
        self.downcast = self.downcast.loc[soak_end:]

    def despike_downcast(self, label, span):
        """
        Still in development. Remove spikes from a data column.
        Adapted from: https://stackoverflow.com/a/71230493

        :param label: string, column label
        :param span: int, moving window size
        :return: Series, de-spiked data column
        """
        fbewma = _ewma_fb(self.downcast[label], span)
        dropped = _drop_outliers(self.downcast[label].to_numpy(), fbewma, deltas[label])
        self.downcast[label] = pd.Series(dropped).interpolate()


def _ewma_fb(column, span):
    """
    Helper function for despike method. Forward/backward exponential weighted
    moving average.

    :param column: array-like, pandas Series
    :param span: int, moving window size
    :return: ndarray
    """
    fwd = pd.Series.ewm(column, span=span).mean()
    bwd = pd.Series.ewm(column[::-1], span=span).mean()
    stacked_ewma = np.vstack((fwd, bwd[::-1]))
    return np.mean(stacked_ewma, axis=0)


def _drop_outliers(column, fbewma, delta):
    """
    Helper function for despike method. Remove outliers from a curve.

    :param column: ndarray, column values
    :param fbewma: ndarray, fwd/bwd weighted moving average
    :param delta: int, outlier limit
    :return: ndarray
    """
    cond_delta = (np.abs(column - fbewma) > delta)
    return np.where(cond_delta, np.nan, column)
