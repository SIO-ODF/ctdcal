import numpy as np
import pandas as pd


def flag_nan_values(data, flag_good=2, flag_nan=9):
    """
    Flag values as either good or nan.

    Parameters
    ----------
    data : array-like
        Variable to be flagged
    flag_good : int, optional
        Flag value for good data
    flag_nan : int, optional
        Flag value for bad (nan) data

    Returns
    -------
    flags : array-like
        Flag for each data point in input
    """

    flags = np.full(np.shape(data), flag_good)
    flags[np.isnan(data)] = flag_nan

    return flags


def flag_outliers(data, flag_good=2, flag_outlier=4, n_sigma=10, ignore_nan=True):
    """
    Flag outliers more than n_sigma standard deviations from the mean.

    Parameters
    ----------
    data : array-like
        Variable to be flagged
    flag_good : int, optional
        Flag value for good data
    flag_outlier : int, optional
        Flag value for outliers
    n_sigma : int, optional
        Number of standard deviations away from mean needed to be outlier
    ignore_nan : bool, optional
        Ignore nan values in data

    Returns
    -------
    flags : array-like
        Flag for each data point in input
    """

    flags = np.full(np.shape(data), flag_good)

    if ignore_nan:
        mean = np.nanmean(data)
        sigma = np.nanstd(data)
    else:
        mean = np.mean(data)
        sigma = np.std(data)

    outliers = np.abs(data - mean) > (n_sigma * sigma)
    flags[outliers] = flag_outlier

    return flags


def stepped_filter(
    residual,
    pressure,
    flag_good=2,
    flag_bad=3,
    threshold=[0.002, 0.005, 0.01, 0.02],
    p_cutoff=[2000, 1000, 500],
):
    """
    Flag residuals using specific thresholds for each pressure range.

    Parameters
    ----------
    residual : array-like
        Residual values
    pressure : array-like
        Pressure values for each data point
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

    residual = np.abs(residual)
    flags = np.full(np.shape(residual), flag_good)

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

    return flags
