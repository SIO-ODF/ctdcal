from ctdcal import flagging
import numpy as np
import pandas as pd
import pytest


# TODO: make a class and recycle 3 data sets, separating outliers and NaNs:
# data = [np.nan] + 97 * [0] + [100, 100]
# data = pd.Series([np.nan] + 97 * [0] + [100, 100])
# data = pd.DataFrame([np.nan] + 97 * [0] + [100, 100])


def test_merge_flags():
    old_flags = 5 * [2]
    new_flags = 5 * [3]

    # check proper flags are kept
    np.testing.assert_array_equal(
        flagging._merge_flags(new_flags, old_flags, keep_higher=True), new_flags
    )
    np.testing.assert_array_equal(
        flagging._merge_flags(new_flags, old_flags, keep_higher=False), old_flags
    )

    # error if flag lengths are different
    with pytest.raises(ValueError):
        flagging._merge_flags(old_flags, new_flags[:-1])
    with pytest.raises(ValueError):
        flagging._merge_flags(old_flags[:-1], new_flags)


def test_nan_values():
    data = [0, 0, 0, np.nan]
    old_flags = len(data) * [2]

    # check NaNs are flagged properly (without old flags)
    assert all(flagging.nan_values(data)[:-1] == 2)
    assert flagging.nan_values(data)[-1] == 9

    # check re-defining flag_good/flag_bad (without old flags)
    new_flags = flagging.nan_values(data, flag_good=6, flag_nan=5)
    assert all(new_flags[:-1] == 6)
    assert new_flags[-1] == 5

    # check NaNs are flagged properly (with old flags)
    assert all(flagging.nan_values(data, old_flags=old_flags)[:-1] == 2)
    assert flagging.nan_values(data, old_flags=old_flags)[-1] == 9

    # check re-defining flag_good/flag_bad (with old flags)
    new_flags = flagging.nan_values(data, old_flags=old_flags, flag_good=6, flag_nan=5)
    assert all(new_flags[:-1] == 6)
    assert new_flags[-1] == 5


def test_outliers():
    data = [np.nan] + 97 * [0] + [100, 100]
    old_flags = (len(data) - 1) * [2] + [7]

    # check outliers are flagged properly (without old flags)
    assert all(flagging.outliers(data)[:-2] == 2)
    assert all(flagging.outliers(data)[-2:] == 4)

    # check re-defining flag_good/flag_bad (without old flags)
    new_flags = flagging.outliers(data, flag_good=6, flag_outlier=5)
    assert all(new_flags[:-2] == 6)
    assert all(new_flags[-2:] == 5)

    # check outliers are flagged properly (with old flags)
    assert all(flagging.outliers(data, old_flags=old_flags)[:-2] == 2)
    assert flagging.outliers(data, old_flags=old_flags)[-1] == old_flags[-1]
    assert flagging.outliers(data, old_flags=old_flags)[-2] == 4

    # check re-defining flag_good/flag_bad (with old_flags)
    new_flags = flagging.outliers(
        data, old_flags=old_flags, flag_good=6, flag_outlier=5
    )
    assert all(new_flags[:-2] == 6)
    assert new_flags[-1] == old_flags[-1]
    assert new_flags[-2] == 5

    # check outliers are not flagged with large n_sigma1
    assert all(flagging.outliers(data, n_sigma1=5)[:-2] == 2)
    assert all(flagging.outliers(data, n_sigma1=5)[-2:] == 4)
    assert all(flagging.outliers(data, n_sigma1=10) == 2)


def test_by_percent_diff():
    data = [np.nan] + 97 * [0] + [100, 100]
    ref = (len(data) - 1) * [0] + [80]
    old_flags = (len(data) - 1) * [2] + [7]

    # check large diffs are flagged properly (without old flags)
    assert all(flagging.by_percent_diff(data, ref)[:-2] == 2)
    assert all(flagging.by_percent_diff(data, ref)[-2:] == 3)

    # check re-defining flag_good/flag_bad (without old flags)
    new_flags = flagging.by_percent_diff(data, ref, flag_good=6, flag_bad=5)
    assert all(new_flags[:-2] == 6)
    assert all(new_flags[-2:] == 5)

    # check large diffs are flagged properly (with old flags)
    new_flags = flagging.by_percent_diff(data, ref, old_flags=old_flags)
    assert all(new_flags[:-2] == 2)
    assert new_flags[-1] == old_flags[-1]
    assert new_flags[-2] == 3

    # check re-defining flag_good/flag_bad (with old flags)
    new_flags = flagging.by_percent_diff(
        data, ref, old_flags=old_flags, flag_good=6, flag_bad=5
    )
    assert all(new_flags[:-2] == 6)
    assert new_flags[-1] == old_flags[-1]
    assert new_flags[-2] == 5

    # check large diffs are not flagged with sufficiently high thresh
    assert flagging.by_percent_diff(data, ref, percent_thresh=10)[-1] == 3
    assert flagging.by_percent_diff(data, ref, percent_thresh=20)[-1] == 2


@pytest.mark.parametrize(
    "offset, p", [(0.003, 2000), (0.006, 1000), (0.011, 500), (0.021, 0)]
)
def test_by_residual(offset, p):
    data = [np.nan] + 97 * [0] + [100, 100]
    ref = np.array(data) + offset
    pres = np.linspace(0, 3000, len(data))
    old_flags = (len(data) - 1) * [2] + [7]

    # check values above threshold are flagged properly (without old flags)
    new_flags = flagging.by_residual(data, ref, pres)
    assert all(new_flags[pres < p] == 2)
    assert all(new_flags[pres > p] == 3)

    # check re-defining flag_good/flag_bad (without old flags)
    new_flags = flagging.by_residual(data, ref, pres, flag_good=6, flag_bad=5)
    assert all(new_flags[pres < p] == 6)
    assert all(new_flags[pres > p] == 5)

    # check values above threshold are flagged properly (with old flags)
    new_flags = flagging.by_residual(data, ref, pres, old_flags=old_flags)
    assert all(new_flags[pres < p] == 2)
    assert all(new_flags[pres > p][:-1] == 3)
    assert new_flags[pres > p][-1] == old_flags[-1]

    # check re-defining flag_good/flag_bad (with old flags)
    new_flags = flagging.by_residual(
        data, ref, pres, old_flags=old_flags, flag_good=6, flag_bad=5
    )
    assert all(new_flags[pres < p] == 6)
    assert all(new_flags[pres > p][:-1] == 5)
    assert new_flags[pres > p][-1] == old_flags[-1]

    # error if threshold not one value longer than p_cutoff
    with pytest.raises(IndexError):
        flagging.by_residual([], [], [], threshold=[1, 2], p_cutoff=[100, 50, 0])
    with pytest.raises(IndexError):
        flagging.by_residual([], [], [], threshold=[1, 2, 3], p_cutoff=[100, 50, 0])
    with pytest.raises(IndexError):
        flagging.by_residual([], [], [], threshold=[1, 2, 3], p_cutoff=[100])

    # error if p_cutoff not in decreasing order
    with pytest.raises(ValueError):
        flagging.by_residual([], [], [], threshold=[1, 2, 3], p_cutoff=[0, 50])
    with pytest.raises(ValueError):
        flagging.by_residual([], [], [], threshold=[1, 2, 3, 4], p_cutoff=[100, 0, 50])
