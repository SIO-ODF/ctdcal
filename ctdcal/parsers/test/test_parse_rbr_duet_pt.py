"""
Unit tests for parse_rbr_duet_pt
"""
from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path
from unittest.mock import patch, Mock

import pandas
import pandas as pd
import pytest

from ctdcal.parsers.parse_rbr_duet_pt import Cast, Parser


class TestParser:
    @pytest.fixture
    def tmp_dir(self, tmp_path):
        path = tmp_path
        fnames = ['a.rsk', 'b.rsk']
        for fname in fnames:
            Path(path, fname).touch()
        return path

class TestCast:
    @pytest.fixture
    def fake_bottle_file(self, tmp_path):
        fake_file = Path(tmp_path, 'spam_btl_mean.pkl')
        pd.DataFrame({'sausage': [0], 'egg': [0], 'spam': [170208000]}).to_pickle(fake_file)
        return fake_file

    def test_load_bottle_times(self, fake_bottle_file):
        cast_id = 'spam'
        cast = Cast(cast_id)

        # test that the parameter is set to something
        cast.load_bottle_times(fake_bottle_file.parent, timecol='spam')
        assert cast.bottle_times is not None
        # test for bad cols
        with pytest.raises(KeyError):
            cast.load_bottle_times(fake_bottle_file.parent, timecol='beans')

    def test_select_raw_file(self):
        cast_id = 'spam'
        cast = Cast(cast_id)
        fake_timestamps = [
                dt.now(),
                dt.now() + td(days=1),
                dt.now() + td(days=2),
                dt.now() + td(days=3)
        ]
        fake_cast_date = fake_timestamps[2] - td(hours=4)

        # test that raw_file is set to the third in the list
        parser = Mock()
        parser.raw_files = ['a', 'm', 'p', 's']
        parser.raw_file_timestamps = fake_timestamps
        parser.start_times = {cast_id: fake_cast_date}
        cast.select_raw_file(parser)
        assert cast.raw_file is 'p'

        # test cast date more recent than raw files
        fake_cast_date = dt.now() + td(days=10)
        parser.start_times = {cast_id: fake_cast_date}
        with pytest.raises(IndexError):
            cast.select_raw_file(parser)

        # test no raw files
        parser.raw_files = []
        with pytest.raises(IndexError):
            cast.select_raw_file(parser)
