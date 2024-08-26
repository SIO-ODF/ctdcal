import json
from pathlib import Path

import pandas as pd
import pytest

from ctdcal.fitting.common import BottleFlags, df_node_to_BottleFlags, get_node, save_node, NodeNotFoundError


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
        # Test non-existing key
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
    def sample_data(self, tmp_path):
        fname = tmp_path / 'data.json'
        with open(fname, 'w') as f:
            json.dump({'spam': {'egg': ['a', 'b'], 'sausage': ['C', 'D']}}, f)
        return fname

    @pytest.fixture
    def sample_df(self):
        spam_dict = {'spam': {'egg': ['a', 'b'], 'sausage': ['C', 'D']}}
        spam_df = pd.DataFrame.from_dict(spam_dict['spam'])
        return spam_df

    def test_df_node_to_BottleFlag(self, sample_df):
        spam = df_node_to_BottleFlags(sample_df)
        assert type(spam.sausage) is list
        assert spam.sausage[0] == 'C'

    def test_get_node(self, sample_data):
        # Test good node
        node = get_node(sample_data, 'spam')
        assert node.sausage[0] == 'C'
        # Test bad node
        with pytest.raises(NodeNotFoundError):
            get_node(sample_data, 'cheese')

    def test_save_node(self, sample_data, tmp_path):
        fname = Path(tmp_path, 'sample_copy.json')
        with open(sample_data, 'r') as f:
            flags = BottleFlags.fromJSON(f.read())
        flags.save(fname)
        spam = get_node(fname, 'spam')
        assert 'x' not in spam.egg
        spam.update_node(egg='x', sausage='Y')
        # Test good node
        save_node(fname, spam, 'spam')
        with open(fname, 'r') as f:
            flags = BottleFlags.fromJSON(f.read())
        assert 'x' in flags.spam.egg
        # Test bad node
        with pytest.raises(NodeNotFoundError):
            save_node(fname, spam, 'limburger')
        save_node(fname, spam, 'Cheddar', create_new=True)
        with open(fname, 'r') as f:
            flags = BottleFlags.fromJSON(f.read())
        assert 'Cheddar' in flags
        # Test no preexisting file
        new_fname = Path(tmp_path, 'sample_new.json')
        assert new_fname.exists() is False
        with pytest.raises(FileNotFoundError):
            save_node(new_fname, spam, 'spam')
        # Test new empty file
        empty_fname = Path(tmp_path, 'sample_empty.json')
        empty_fname.touch()
        save_node(empty_fname, spam, 'spam', create_new=True)
        with open(empty_fname, 'r') as f:
            flags = BottleFlags.fromJSON(f.read())
        assert 'spam' in flags
