"""
Classes, definitions and utilities for all ctdcal modules
"""
import logging
from io import BufferedIOBase, BytesIO, StringIO
from pathlib import Path
from typing import Union
from zipfile import is_zipfile, ZipFile
from zipimport import ZipImportError

import pandas as pd
import requests
import yaml

from munch import munchify

from ctdcal import get_ctdcal_config

log = logging.getLogger(__name__)
cfg = get_ctdcal_config()


# Constants
# ---------


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
    cfgfile = validate_file(cfgfile)
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


# File imports
def load_cnv(cnv_file: Union[str, Path]) -> pd.DataFrame:
    """
    Load Sea-Bird converted (.cnv) cast file into DataFrame
    TODO: Union is no longer required in python 3.10. use `cnv_file: str | Path` instead.
    """
    with open(cnv_file) as f:
        file = f.readlines()

    # parse column names
    info = dict()
    cols = []
    for idx, line in enumerate(file):
        # get variable info
        if line.strip("# \n").startswith(("nquan", "nvalues", "units", "bad_flag")):
            k, v = line.strip("# \n").split("=")
            info[k.strip()] = v.strip()

        # get column names
        elif line.strip("# ").startswith("name"):
            # expected format is:   # name 0 = col_name: long_description
            cols.append(line.split(":")[0].split("=")[-1].strip())

        # last row before data begins
        elif line.startswith("*END*"):
            data_index = idx + 1
            break

        # anything else is a comment line
        else:
            continue  # pragma: no cover

    # read data
    return pd.read_csv(
        cnv_file,
        skiprows=range(0, data_index),
        # delim_whitespace=True,
        sep=r'\s+',
        names=cols,
        engine="python",
        skipinitialspace=True,
        na_values=info["bad_flag"],
    )


def load_exchange_btl(btl_file: Union[str, Path]) -> pd.DataFrame:
    """
    Load WHP-exchange bottle file (_hy1.csv) into DataFrame. File can be on local
    file system or downloaded from an appropriate cchdo.ucsd.edu link
    (e.g., https://cchdo.ucsd.edu/data/19436/325020210316_hy1.csv)

    Adapted from cchdo.hydro package.

    TODO: Union is no longer required in python 3.10. use `|` instead.

    Parameters
    ----------
    btl_file : str or Path
        Name or URL of file to be loaded

    Returns
    -------
    df : DataFrame
        Loaded bottle file
    """
    # read from url
    if isinstance(btl_file, (str, Path)) and str(btl_file).startswith("http"):
        log.info(f"Loading bottle file {Path(btl_file).name} from http link")
        file = requests.get(btl_file).text.splitlines(keepends=True)

    # read from file
    elif isinstance(btl_file, (str, Path)):
        log.info(f"Loading bottle file {Path(btl_file).name} from local file")
        with open(btl_file) as f:
            file = f.readlines()

    # find index of units row
    for idx, line in enumerate(file):
        # skip comment lines (which may reference EXPOCODE and break membership test)
        if line.strip().startswith("#"):
            continue

        # find index of units row
        if "EXPOCODE" in line:
            units = idx + 1  # units row immediately follows column names
            break

    return pd.read_csv(
        StringIO("".join(file)),
        skiprows=[0, units],
        skipfooter=1,
        engine="python",
        comment="#",
        skipinitialspace=True,
    )


# File exports
def list_to_file(fname, outdir, lst):
    """
    Writes a list to a file, one list item to a line.

    Parameters
    ----------
    fname : str or Path-like
    outdir : str or Path-like
    lst : list
    """
    outdir = validate_dir(outdir)
    outfile = Path(outdir, fname)
    with open(outfile, 'w') as f:
        f.write('\n'.join(lst))

def load_exchange_ctd(
    ctd_file: Union[str, Path, BufferedIOBase],
    n_files=None,
    recursed=False,
) -> pd.DataFrame:
    """
    Load WHP-exchange CTD file(s) (_ct1.csv) into DataFrame. File(s) can be on local
    file system or downloaded from an appropriate cchdo.ucsd.edu link
    (e.g., https://cchdo.ucsd.edu/data/19434/325020210316_ct1.zip)

    Adapted from cchdo.hydro package.

    Parameters
    ----------
    ctd_file : str or Path
        Name or URL of file to be loaded

    n_files : int, optional
        Number of files to load from .zip archive

    Returns
    -------
    header : dict or list of dict
        File metadata from header(s) (e.g., EXPOCODE, STNNBR, CASTNO)
    df : DataFrame or list of DataFrame
        Loaded CTD file(s)
    TODO: Union is no longer required in python 3.10. use `|` instead.
    """
    # read from url (.zip)
    if isinstance(ctd_file, (str, Path)) and str(ctd_file).startswith("http"):
        log.info(f"Loading CTD file {Path(ctd_file).name} from http link")
        data_raw = BytesIO(requests.get(ctd_file).content)

    # read from file
    elif isinstance(ctd_file, (str, Path)):
        log.info(f"Loading CTD file {Path(ctd_file).name} from local file")
        with open(ctd_file, "rb") as f:
            data_raw = BytesIO(f.read())

    # read from open file
    elif isinstance(ctd_file, BufferedIOBase):
        log.info("Loading open file object")
        data_raw = BytesIO(ctd_file.read())

    # .zip special behavior
    if is_zipfile(data_raw):
        log.info("Loading CTD files from .zip")

        if recursed is True:
            raise ZipImportError("Recursive .zip files encountered... exiting")

        data_raw.seek(0)  # is_zipfile moves cursor to EOF, reset to start
        zip_contents = []
        with ZipFile(data_raw) as zf:
            for zipinfo in zf.infolist():
                zip_contents.append(BytesIO(zf.read(zipinfo)))

        # list comprehension is same as using functools.partial, just different syntax
        return zip(
            *[load_exchange_ctd(zc, recursed=True) for zc in zip_contents[:n_files]]
        )

    else:
        data_raw.seek(0)  # is_zipfile moves cursor to EOF, reset to start
        file = data_raw.read().decode("utf8").splitlines(keepends=True)

    # process metadata
    for idx, line in enumerate(file):
        # skip comment lines (which may reference CTDPRS and break membership test)
        if line.strip().startswith("#"):
            continue

        # find header info
        if line.startswith("NUMBER_HEADERS"):
            header_ind = idx

        # find index of units row
        if "CTDPRS" in line:
            columns = idx
            units = idx + 1  # units row immediately follows column names
            break

    # break down header rows
    header = {}
    for line in file[header_ind:columns]:
        k, v = line.strip("\n").split("=")
        header[k.strip()] = v.strip()

    return header, pd.read_csv(
        StringIO("".join(file)),
        skiprows=list(range(0, columns)) + [units],  # skip up to column names (+ units)
        skipfooter=1,
        engine="python",
        comment="#",
        skipinitialspace=True,
    )


## Cast list functions
def make_cast_id_list(rawdir, outdir, fname="cast_id_list.csv"):
    """
    Attempt to automatically generate list of station/casts from raw files.
    """
    search_dir = validate_dir(Path(rawdir))
    raw_files = Path(search_dir).glob("*.hex")
    cast_id_list = sorted([f.stem for f in raw_files])
    if len(cast_id_list) < 1:
        raise FileNotFoundError('No raw data files found.')
    # pd.Series(ssscc_list, dtype=str).to_csv(fname, header=None, index=False, mode="x")
    outdir = validate_dir(outdir, create=True)
    list_to_file(fname, outdir, cast_id_list)
    return cast_id_list



def get_cast_id_list(fname, rawdir, outdir, auto_generate=True):
    """
    Loads a cast list from a file. If the file is not found, the cast list will be
    auto-generated from a directory of raw files when auto_generate is True (default
    behavior).

    Parameters
    ----------
    fname
    auto_generate

    Returns
    -------

    """
    fname = Path(outdir, fname)
    try:
        # test if file exists
        fname = validate_file(fname)
    except FileNotFoundError:
        # auto-generate or die
        if auto_generate is True:
            make_cast_id_list(rawdir, outdir, fname=fname)
        else:
            raise

    # file exists or has been auto-generated
    with open(fname, "r") as f:
        # skip comment lines
        id_list = [line.strip() for line in f.readlines() if not line.startswith('#')]
    return id_list

def get_ssscc_list(fname="data/ssscc.csv"):
    """
    Load a list of casts from a file.

    Parameters
    ----------
    fname : path_like
        Input file. Type is anything that can be interpreted by Python as a
        path, such as a string or a Pathlib object.

    Returns
    -------
    list
        Cast names or identifiers, as a list of strings.
    """
    log.warning("Use of get_ssscc_list() is deprecated. Use get_cast_id_list() instead.")
    ssscc_list = []
    with open(fname, "r") as lines:
        for line in lines:
            # skip comment lines
            if not line.startswith("#"):
                ssscc_list.append(line.strip())
    return ssscc_list


# Utilities

def print_progress_bar(
        iteration,
        total,
        prefix="",
        suffix="",
        decimals=1,
        length=100,
        fill="█",
        printEnd="\r",
):
    """
    A progress bar, helpful for implementing into loops or highlighting progression through processing.

    https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters/13685020
    credit: u/Greenstick
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)  # Potential to add to log
    # Print New Line on Complete
    if iteration == total:
        print()
