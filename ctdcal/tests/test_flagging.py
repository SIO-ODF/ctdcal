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
    assert all(flagging.nan_values(data, flag_good=6)[:-1] == 6)
    assert flagging.nan_values(data)[-1] == 9
    assert flagging.nan_values(data, flag_nan=5)[-1] == 5

    # check NaNs are flagged properly (with old flags)
    assert all(flagging.nan_values(data, old_flags=old_flags)[:-1] == 2)
    assert all(flagging.nan_values(data, old_flags=old_flags, flag_good=6)[:-1] == 6)
    assert flagging.nan_values(data, old_flags=old_flags)[-1] == 9
    assert flagging.nan_values(data, old_flags=old_flags, flag_nan=5)[-1] == 5


def test_outliers():
    data = [np.nan] + 97 * [0] + [100, 100]
    # old_flags = len(data) * [2]

    assert all(flagging.outliers(data)[:-2] == 2)
    assert all(flagging.outliers(data)[-2:] == 4)
