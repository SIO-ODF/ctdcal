from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

import ctdcal.common
from ctdcal.common import (
    load_user_config,
    validate_dir,
    validate_file,
    make_cast_id_list,
    get_cast_id_list,
    load_fit_groups,
)


class TestUserConfig:
    """
    Tests and fixtures for working with user configurations.
    """
    @pytest.fixture
    def tmp_cfg(self, tmp_path):
        fname = tmp_path / 'cfg.yaml'
        with open(fname, 'w') as f:
            yaml.dump({'cheese': {'cheddar': 0, 'limburger': 1}}, f)
        return fname

    def test_load_user_config(self, tmp_cfg):
        # Test good cfgfile
        cfg = load_user_config(tmp_cfg)
        assert cfg.cheese.limburger == 1

        # Test bad cfgfile
        with pytest.raises(FileNotFoundError):
            load_user_config('fake_file')


class TestPathValidation:
    """
    Unit tests and fixtures for file and directory-related validation
    functions.
    """
    @pytest.fixture
    def tmp_dir(self, tmp_path):
        path = tmp_path
        return path

    def test_validate_dir(self, tmp_dir):
        p = 'spam/'
        # Test for bad input
        with pytest.raises(TypeError):
            validate_dir(None)

        # Test for bad path
        with pytest.raises(FileNotFoundError):
            validate_dir(p, create=False)

        # Test for good path
        with patch('pathlib.Path.is_dir', return_value=True):
            validated_path = validate_dir(p, create=False)
        assert validated_path == Path(p)

        # Test that bad dir is created
        test_dir = tmp_dir / p / 'egg'
        assert test_dir.is_dir() is False
        assert validate_dir(test_dir, create=True).is_dir()

        # Test for path exists but is not a dir
        samename = tmp_dir / p / 'cheese'
        assert samename.exists() is False
        samename.touch()
        with pytest.raises(FileExistsError):
            validate_dir(samename, create=False)
        with pytest.raises(FileExistsError):
            validate_dir(samename, create=True)

    def test_validate_file(self, tmp_dir):
        p = 'spam.file'
        # Test for bad input
        with pytest.raises(TypeError):
            validate_file(None)

        # Test for bad path
        with pytest.raises(FileNotFoundError):
            validate_file(p, create=False)

        # Test for good path
        with patch('pathlib.Path.is_file', return_value=True):
            validated_path = validate_file(p, create=False)
        assert validated_path == Path(p)

        # Test that non-existent file is created
        tmp_file = tmp_dir / 'spam' / 'egg.file'
        assert tmp_file.is_file() is False
        assert validate_file(tmp_file, create=True).is_file()

        # Test that filename exists but is not a file
        samename = tmp_dir / 'spam' / 'egg'
        assert samename.exists() is False
        samename.mkdir(parents=True)
        with pytest.raises(FileExistsError):
            validate_file(samename, create=False)
        with pytest.raises(FileExistsError):
            validate_file(samename, create=True)


class TestCastListFunctions:
    """
    Unit tests and fixtures for functions which generate or manipulate
    cast file lists.
    """
    @pytest.fixture
    def tmp_file(self, tmp_path):
        fname = tmp_path / 'casts.csv'
        cast_list = ['a', 'b', 'c', 'd', 'f']
        with open(fname, 'w') as f:
            f.write('\n'.join(cast_list))
        return fname

    @pytest.fixture
    def tmp_dir(self, tmp_path):
        path = tmp_path
        return path

    def test_make_cast_id_list(self, tmp_dir):
        inst = 'inst'
        validate_dir(Path(tmp_dir, 'raw', inst), create=True)
        # test happy path
        with patch('pathlib.Path.glob', return_value=[Path(p) for p in ['a.spam', 'b.spam', 'f.spam']]):
            id_list = make_cast_id_list(tmp_dir, inst, tmp_dir)
        assert len(id_list) == 3
        assert id_list[2] == 'f'

        # test no raw files
        with patch('pathlib.Path.glob', return_value=[]):
            with pytest.raises(FileNotFoundError):
                make_cast_id_list(tmp_dir, inst, tmp_dir)

    @patch('ctdcal.common.make_cast_id_list')
    def test_get_cast_id_list(self, mock_make, tmp_dir, tmp_file):
        # test happy path
        id_list = get_cast_id_list(tmp_file, None, None, tmp_dir)
        assert len(id_list) == 5
        assert id_list[2] == 'c'

        # test cast id file not found
        with pytest.raises(FileNotFoundError):
            get_cast_id_list('fake_file', None, None, tmp_dir, auto_generate=False)

        # test auto-generate
        def side_effect(outdir):
            fake_file = outdir / 'fake_file'
            ctdcal.common.list_to_file(fake_file, tmp_dir, ['a', 'b', 'f'])
            return None

        mock_make.side_effect = side_effect(tmp_dir)
        id_list = get_cast_id_list('fake_file', None, None, tmp_dir)
        assert len(id_list) == 3
        assert id_list[2] == 'f'

    def test_load_fit_groups(self):
        fake_groups = {'spam': ['a', 'd']}
        fake_list = ['a', 'b', 'c', 'd', 'f']

        # happy path
        # spam has two groups:
        fit_groups_all = load_fit_groups(fake_list, fake_groups)
        assert len(fit_groups_all) == 1
        assert len(fit_groups_all['spam']) == 2
        assert len(fit_groups_all['spam'][0]) == 3
        assert 'd' not in fit_groups_all['spam'][0]
        assert 'd' in fit_groups_all['spam'][1]

        # spam is all one group:
        fake_groups = {'spam': ['a']}
        fit_groups_all = load_fit_groups(fake_list, fake_groups)
        assert len(fit_groups_all['spam']) == 1
        assert len(fit_groups_all['spam'][0]) == 5

        # sad path
        # ambiguous slice indexer (duplicate items in list)
        fake_groups = {'spam': ['a', 'd']}
        fake_list = ['a', 'b', 'd', 'd', 'f']
        with pytest.raises(ValueError):
            fit_groups_all = load_fit_groups(fake_list, fake_groups)

        # slice indexer not found
        fake_list = ['a', 'b', 'c', 'f']
        with pytest.raises(ValueError):
            fit_groups_all = load_fit_groups(fake_list, fake_groups)
