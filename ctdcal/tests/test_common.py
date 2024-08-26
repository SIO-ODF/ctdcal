from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from ctdcal.common import load_user_config, validate_dir, validate_file


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
        # TODO: Test for invalid file?


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
