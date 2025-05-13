"""
Classes and functions for flagging data by all modules.
"""
import logging

import numpy as np
import pandas as pd

from ctdcal import get_ctdcal_config
from ctdcal.plotting.plot_fit import _intermediate_residual_plot

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def _flag_btl_data(
    df,
    param=None,
    ref=None,
    cast_id='SSSCC',
    thresh=[0.002, 0.005, 0.010, 0.020],
    f_out=None,
):
    """
    Flag CTD "btl" data against reference measurement (e.g. SBE35, bottle salts).

    Parameters
    ----------
    df : DataFrame,
        DataFrame containing btl data
    param : str
        Name of parameter to calibrate (e.g. "CTDCOND1", "CTDTMP2")
    ref : str
        Name of reference parameter to calibrate against (e.g. "BTLCOND", "T90")
    thresh : list of float, optional
        Maximum acceptable residual for each pressure range
    f_out : str, optional
        Path and filename to save residual vs. pressure plots

    Returns
    -------
    df_ques : DataFrame
        Data flagged as questionable (flag 3s)
    df_bad : DataFrame
        Data flagged as bad (flag 4s)

    """
    # prs = cfg.column["p"]

    # Remove extreme outliers and code bad
    df = df.reset_index(drop=True)
    df["Diff"] = df[ref] - df[param]
    df["Flag"] = outliers(df["Diff"])

    # Find values that are above the threshold and code questionable
    df["Flag"] = by_residual(df[param], df[ref], df['CTDPRS'], old_flags=df["Flag"])
    df_good = df[df["Flag"] == 2].copy()
    df_ques = df[df["Flag"] == 3].copy()
    df_bad = df[df["Flag"] == 4].copy()

    if f_out is not None:
        if param == cfg.column["t1"]:
            xlabel = "T1 Residual (T90 C)"
        elif param == cfg.column["t2"]:
            xlabel = "T2 Residual (T90 C)"
        elif param == cfg.column["c1"]:
            xlabel = "C1 Residual (mS/cm)"
        elif param == cfg.column["c2"]:
            xlabel = "C2 Residual (mS/cm)"
        f_out = str(f_out).split(".pdf")[0] + "_postfit.pdf"
        _intermediate_residual_plot(
            df["Diff"],
            df['CTDPRS'],
            df[cast_id],
            show_thresh=True,
            xlabel=xlabel,
            f_out=f_out,
        )
        f_out = str(f_out).split(".pdf")[0] + "_flag2.pdf"
        _intermediate_residual_plot(
            df_good["Diff"],
            df_good['CTDPRS'],
            df_good[cast_id],
            show_thresh=True,
            xlabel=xlabel,
            f_out=f_out,
        )

    return df_ques, df_bad


def _merge_flags(new_flags, old_flags, keep_higher=True):
    """
    Merge old and new flags into one flag vector.

    Returns new_flags if old_flags is None.

    Parameters
    ----------
    new_flags : array-like
        Newly calculated flags
    old_flags : int, optional
        User-supplied original data flags
    keep_higher : bool, optional
        Whether to keep higher (default) or lower value flags

    Returns
    -------
    merged_flags : array-like
        Merged old and new flags
    """

    if old_flags is None:
        return new_flags

    new_flags = np.squeeze(new_flags)
    old_flags = np.squeeze(old_flags)
    merged_flags = np.copy(new_flags)

    if keep_higher:
        is_higher = old_flags > new_flags
        merged_flags[is_higher] = old_flags[is_higher]
    else:
        # can't really think of a use case for this... maybe delete?
        is_lower = old_flags < new_flags
        merged_flags[is_lower] = old_flags[is_lower]

    return merged_flags


def nan_values(data, old_flags=None, flag_good=2, flag_nan=9):
    """
    Flag values as either good or nan.

    Parameters
    ----------
    data : array-like
        Variable to be flagged
    old_flags : array-like, optional
        Original data flags to be merged in (if provided)
    flag_good : int, optional
        Flag value for good data
    flag_nan : int, optional
        Flag value for bad (nan) data

    Returns
    -------
    flags : array-like
        Flag for each data point in input
    """

    data = np.squeeze(data)
    flags = np.full(np.shape(data), flag_good)
    flags[np.isnan(data)] = flag_nan

    return _merge_flags(flags, old_flags)


def outliers(
    data,
    old_flags=None,
    flag_good=2,
    flag_outlier=4,
    n_sigma1=2,
    n_sigma2=20,
    ignore_nan=True,
):
    """
    Flag extreme outliers using standard deviations from the mean as a threshold.

    Outliers are identified over two passes. For the first pass, mean and standard
    deviation of data are calculated for all data. Values more than n_sigma1 standard
    deviations from mean are (temporarily) flagged questionable. For the second pass,
    mean and standard deviation are re-calculated with questionable data excluded. Data
    more than n_sigma2 standard deviations from mean are flagged as outliers.

    Parameters
    ----------
    data : array-like
        Variable to be flagged
    old_flags : array-like, optional
        Original data flags to be merged in (if provided)
    flag_good : int, optional
        Flag value for good data
    flag_outlier : int, optional
        Flag value for outliers
    n_sigma1 : int, optional
        Number of standard deviations away from mean needed to be excluded from statistics
    n_sigma2 : int, optional
        Number of standard deviations away from mean needed to be outlier
    ignore_nan : bool, optional
        Ignore nan values in data

    Returns
    -------
    flags : array-like
        Flag for each data point in input

    Notes
    -----
    Functionality is similar to Sea-Bird's "Wild Edit" in Seasoft V2.
    """

    data = np.squeeze(data)

    # np.squeeze converts DataFrame->Series but does nothing to Series
    if type(data) is pd.Series:
        data = data.to_numpy()  # force to np.array

    flags = np.full(np.shape(data), flag_good).squeeze()

    # function aliases
    if ignore_nan:
        mean, std = np.nanmean, np.nanstd
    else:
        mean, std = np.mean, np.std

    # pass 1
    questionable = np.abs(data - mean(data)) > (n_sigma1 * std(data))

    # pass 2
    data_mean = mean(data[~questionable])
    data_std = std(data[~questionable])
    outliers = np.abs(data - data_mean) > (n_sigma2 * data_std)
    flags[outliers] = flag_outlier

    return _merge_flags(flags, old_flags)


def by_percent_diff(
    data, ref_data, old_flags=None, percent_thresh=1, flag_good=2, flag_bad=3
):
    """
    Flag residuals based on the percent error relative to the reference data.

    Parameters
    ----------
    data : array-like
        Variable to be flagged
    ref_data : array-like
        Reference data to compare against
    old_flags : array-like, optional
        Original data flags to be merged in (if provided)
    percent_thresh : int, optional
        Percent different cutoff to be flagged as bad data
    flag_good : int, optional
        Flag value for good data
    flag_bad : int, optional
        Flag value for bad data

    Returns
    -------
    flags : array-like
        Flag for each data point in input
    """
    data = np.squeeze(data)
    ref_data = np.squeeze(ref_data)

    #   Avoid dividing by zero whenever possible
    data = np.where(data == 0, np.nan, data)

    percent_diff = (np.abs(data - ref_data) / data).squeeze() * 100
    flags = np.full(np.shape(percent_diff), flag_good).squeeze()

    flags[percent_diff > percent_thresh] = flag_bad

    return _merge_flags(flags, old_flags)


def by_residual(
    data,
    ref_data,
    pressure,
    old_flags=None,
    flag_good=2,
    flag_bad=3,
    threshold=[0.002, 0.005, 0.01, 0.02],
    p_cutoff=[2000, 1000, 500],
):
    """
    Flag residuals using specific thresholds for each pressure range.

    Parameters
    ----------
    data : array-like
        Variable to be flagged
    ref_data : array-like
        Reference data to compare against
    flag_good : int, optional
        Flag value for good data
    flag_bad : int, optional
        Flag value for data outside threshold
    threshold : list of float, optional
        Threshold between good and bad values
    p_cutoff : list of int, optional
        Edges of pressure bins

    Returns
    -------
    flags : array-like
        Flag for each data point in input
    """
    if (len(threshold) - len(p_cutoff)) != 1:
        raise IndexError("threshold must be one value longer than p_cutoff")

    if any(np.diff(p_cutoff) >= 0):
        raise ValueError("p_cutoff must be monotonically decreasing")

    residual = np.abs(data - ref_data).squeeze()
    flags = np.full(np.shape(residual), flag_good).squeeze()

    # first and last pressure bin need special handling
    is_deep = pressure > p_cutoff[0]
    is_shallow = pressure < p_cutoff[-1]
    flags[is_deep & (residual > threshold[0])] = flag_bad
    flags[is_shallow & (residual > threshold[-1])] = flag_bad

    p_bins = zip(p_cutoff, p_cutoff[1:])  # make pairs (e.g. (2000, 1000), (1000, 500))
    threshold = threshold[1:-1]  # deep and shallow threshold values
    for (p_upper, p_lower), thresh in zip(p_bins, threshold):
        in_bin = (pressure <= p_upper) & (pressure > p_lower)
        flags[in_bin & (residual > thresh)] = flag_bad

    return _merge_flags(flags, old_flags)
