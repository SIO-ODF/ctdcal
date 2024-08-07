import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from ctdcal.fitting.common import BottleFlags, df_node_to_BottleFlag


class TestBottleFlags:
    @pytest.fixture
    def sample_data(self, tmp_path):
        fname = tmp_path / 'data.json'
        with open(fname, 'w') as f:
            json.dump({'spam': {'egg': []}}, f)
        return fname

    @pytest.fixture
    def sample_fname(self, tmp_path):
        fname = tmp_path / 'save.json'
        return fname

    def test_add_node(self):
        flags = BottleFlags()
        flags.add_node('spam', ['egg', 'baked beans'])
        assert 'spam' in flags
        assert 'baked beans' in flags.spam
        assert len(flags.spam.egg) == 0

    def test_update_node(self, sample_data):
        with open(sample_data, 'r') as f:
            flags = BottleFlags.fromJSON(f.read())
        assert 'lovely spam' not in flags.spam.egg
        # Test existing key
        flags.spam.update_node(egg='lovely spam')
        assert 'lovely spam' in flags.spam.egg
        # Test nonexisting key
        with pytest.raises(KeyError):
            flags.spam.update_node(cheese='cheddar')

    def test_save(self, sample_data, sample_fname):
        with open(sample_data, 'r') as f:
            flags = BottleFlags.fromJSON(f.read())
        assert not sample_fname.is_file()
        flags.save(sample_fname)
        assert sample_fname.is_file()


class TestBottleFlagWrangling:
    @pytest.fixture
    def sample_df(self):
        spam_dict = {'spam': {'egg': ['a', 'b'], 'sausage': ['C', 'D']}}
        spam_df = pd.DataFrame.from_dict(spam_dict['spam'])
        return spam_df

    def test_df_node_to_BottleFlag(self, sample_df):
        bf = df_node_to_BottleFlag(sample_df, 'spam')
        assert type(bf.spam.sausage) is list
        assert bf.spam.sausage[0] == 'C'
