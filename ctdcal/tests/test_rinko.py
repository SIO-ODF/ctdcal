import numpy as np

from ctdcal import rinko


def test_salinity_correction():

    DO_c = np.array([100, 100, 0, 0, 0, np.nan])
    T = np.array([1, 1, 0, 0, np.nan, 0])
    C = np.array([1, 1, 0, np.nan, 0, 0])

    # check values are converted correctly
    corr = rinko.salinity_correction(DO_c, T, C)
    assert corr.shape == (6,)
    assert corr[0] == corr[1]
    assert corr[2] == 0
    assert all(np.isnan(corr[3:]))
