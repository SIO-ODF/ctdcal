"""
Classes, definitions and utilities for all ctdcal modules
"""
from pathlib import Path
import yaml

from munch import munchify


# Function definitions
# --------------------

# Configuration
def load_user_config(cfgfile):
    """
    Load user-defined parameters from a configuration file. Return a Munch
    object (dictionary).

    Parameters
    ----------
    cfgfile : str or Path-like
        Path to the configuration file.

    Returns
    -------
    Munch object
    """
    with open(cfgfile, 'r') as f:
        cfg = yaml.safe_load(f)
        return munchify(cfg)


# Input Validation
def validate_dir(pathname, create=False):
    """
    Test if a directory exists, and optionally create it if it does not. Raise
    an exception if the directory does not exist, cannot be created, or is not
    a directory. Return the validated path.

    Parameters
    ----------
    pathname : str or PathLike
        Directory to validate.
    create : bool
        If true, create the directory if it does not exist. Default is false.

    Returns
    -------
    Path object
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
    Test if a file exists, and optionally create it if it does not. Raise
    an exception if the file does not exist, cannot be created, or is not
    a file. Return the validated path.

    Parameters
    ----------
    pathname : str or Path-like
        Filename to validate.
    create : bool
        If true, create the file if it does not exist. Default is false.

    Returns
    -------
    Path object
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
