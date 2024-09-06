"""
Unit tests for ctdcal.processors.cast_tools
"""
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from ctdcal.processors.cast_tools import Cast


class TestCast:
    @pytest.fixture
    def tmp_df(self, tmp_path):
        data = {'spam': [1, 2, 3, 2, 1, 2, 3, 4, 5, 4, 3, 2, 1]}
        return pd.DataFrame(data=data)

    @pytest.fixture
    def rnd_df(self):
        cols = ['scan_datetime', 'P', 'GPSLAT', 'GPSLON', 'ALT']
        data = np.random.randint(10, size=[4, 5])
        return pd.DataFrame(data, columns=cols)

    def test_load_cast(self, tmp_df):
        # test happy path
        with patch('pandas.read_pickle', return_value=tmp_df):
            cast = Cast('cast', 'fake_path')
        assert type(cast.proc) == pd.DataFrame
        # test sad path
        with pytest.raises(FileNotFoundError):
            Cast('spam', 'fake_path')

    def test_parse_downcast(self, tmp_df):
        with patch('pandas.read_pickle', return_value=tmp_df):
            cast = Cast('cast', 'fake_dir')
        # test p_col not set
        with pytest.raises(AttributeError):
            cast.parse_downcast('fake_data')
        # test happy path
        cast.p_col = 'spam'
        assert cast.downcast is None
        cast.parse_downcast(cast.proc)
        assert type(cast.downcast) == pd.DataFrame
        assert cast.downcast is not None
        assert cast.downcast['spam'].values[-1] == 5

    def test_filter(self, tmp_df):
        with patch('pandas.read_pickle', return_value=tmp_df):
            cast = Cast('spam', 'fake_dir')
        cast.filter(cast.proc, 2, cols=['spam'])
        assert type(cast.filtered) == pd.DataFrame
        # test bad filter type
        with pytest.raises(AttributeError):
            cast.filter(cast.proc, 2, win_type='egg')

    def test_trim_soak(self, tmp_df):
        with patch('pandas.read_pickle', return_value=tmp_df):
            cast = Cast('spam', 'fake_dir')
        # test p_col not set
        with pytest.raises(AttributeError):
            cast.trim_soak('fake_data', 'fake_win')
        # test downcast not parsed
        cast.p_col = 'spam'
        with pytest.raises(AttributeError):
            cast.trim_soak('fake_data', 'fake_win')
        # test happy path
        cast.parse_downcast(cast.proc)
        cast.trim_soak(cast.downcast, 3, max_soak=2)
        assert type(cast.trimmed) == pd.DataFrame

    def test_get_details(self, rnd_df):
        with patch('pandas.read_pickle', return_value=rnd_df):
            cast = Cast('spam', 'fake_dir')
            # test p_col not set
            with pytest.raises(AttributeError):
                cast.get_details()
        # test downcast not parsed and trimmed
        cast.p_col = 'P'
        with pytest.raises(AttributeError):
            cast.get_details()
        # test happy path
        cast.trimmed = cast.proc
        spam = cast.get_details()
        assert spam.size == 9
        # test no altimeter
        cast.proc.drop('ALT', axis=1)
        cast.trimmed = cast.proc
        spam = cast.get_details()
        assert 'ALT' not in spam.columns
        # test incomplete gps
        cast.proc.drop('GPSLAT', axis=1)
        spam = cast.get_details()
        assert 'GPSLON' not in spam.columns

    # def test_get_pressure_offsets(self):
    #     # TODO: write this when get_pressure_offsets is updated
    #     assert False
