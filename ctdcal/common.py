#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:package: ctdcal.common
:file: ctdcal/common.py
:author: Allen Smith
:brief: Classes, definitions and utilities for all ctdcal modules
"""
from pathlib import Path
import yaml

from munch import munchify


# Function definitions
# --------------------

# Configuration
def load_user_config(cfgfile):
    with open(cfgfile, 'r') as f:
        cfg = yaml.safe_load(f)
        return munchify(cfg)

# Input Validation
def validate_dir(pathname, create=False):
    """
    Test if a directory exists, and optionally create it if it does not. Raises
    an exception if the directory does not exist, cannot be created, or is not
    a directory.

    Parameters
    ----------
    pathname - (str, PathLike) directory to validate
    create - (bool) create the directory if true

    Returns
    -------
    Path object for the validated directory
    """
    p = Path(pathname)
    if create is True:
        p.mkdir(parents=True, exist_ok=True)
        return p
    elif p.is_dir():
        return p
    elif p.exists():
        raise FileExistsError("%s already exists but is not a directory." % str(p))
    else:
        raise FileNotFoundError("The directory %s could not be found" % str(p))


def validate_file(pathname, create=False):
    """
    Test if a file exists, and optionally create it if it does not. Raises
    an exception if the file does not exist, cannot be created, or is not
    a file.

    Parameters
    ----------
    pathname - (str, PathLike) filename to validate
    create - (bool) create the file if true

    Returns
    -------
    Path object for the validated file
    """
    p = Path(pathname)
    if create is True:
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.is_dir():
            raise FileExistsError
        p.touch(exist_ok=True)
        return p
    elif p.is_file():
        return p
    elif p.exists():
        raise FileExistsError("%s exists but is not a file." % str(p))
    else:
        raise FileNotFoundError("The file %s could not be found." % str(p))
