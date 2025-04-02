"""
Processes reference temperature data.
"""
import csv
import logging
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

from ctdcal import get_ctdcal_config
from ctdcal.common import validate_dir
from ctdcal.parsers.parse_sbe35 import parse_sbe35

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def reft_loader(infile):
    cols = ['index_memory', 'datetime', 'btl_fire_num', 'diff', 'raw_value', 'T90']
    dtypes = {
            'index_memory': int,
            'datetime': object,
            'btl_fire_num': int,
            'diff': int,
            'raw_value': float,
            'T90': float
    }
    # read csv into a dataframe
    reft_df = pd.read_csv(infile, header=0, names=cols, dtype=dtypes)
    # convert to native datetime format
    reft_df['datetime'] = pd.to_datetime(reft_df['datetime'])

    # assign initial flags (large "diff" = unstable reading, flag questionable)
    reft_df["REFTMP_FLAG_W"] = 2
    reft_df.loc[reft_df["diff"].abs() >= 3000, "REFTMP_FLAG_W"] = 3

    return reft_df

def _reft_loader(ssscc, reft_dir=cfg.dirs["reft"]):
    """
    Loads SBE35.cap files and assembles into a dataframe.

    Parameters
    ----------
    ssscc : str
        Station to load .CAP file for
    reft_dir : str, pathlib Path
        Path to the reft folder

    Returns
    -------
    reftDF : DataFrame
        DataFrame of .CAP file with headers

    """
    log.warning("Use of _reft_loader() is deprecated. Use reft_loader() instead.")
    # semi-flexible search for reft file (in the form of *ssscc.cap)
    try:
        reft_path = sorted(Path(reft_dir).glob(f"*{ssscc}.cap"))[0]
    except IndexError:
        raise FileNotFoundError

    # this works better than pd.read_csv as format is semi-inconsistent (cf .cap files)
    with open(reft_path, "r", newline="") as f:
        reftF = csv.reader(
            f, delimiter=" ", quoting=csv.QUOTE_NONE, skipinitialspace="True"
        )
        reftArray = []
        for row in reftF:
            if len(row) != 17:  # skip over 'bad' rows (empty lines, comments, etc.)
                if len(row) > 17:
                    log.warning(f"Raw REFT for {ssscc} includes line with abnormally large number of columns. Check .CAP file.")
                else:
                    continue
            reftArray.append(row)
    if len(reftArray)>36:
        log.warning("Raw REFT file for {ssscc} exceeds 36 entries. Check .CAP file.")

    reftDF = pd.DataFrame.from_records(reftArray)
    pd.set_option('future.no_silent_downcasting', True) #   Opt in to future downcasting
    reftDF = reftDF.replace(
        to_replace=["bn", "diff", "val", "t90", "="], value=np.nan
    )
    reftDF = reftDF.dropna(axis=1)
    reftDF.loc[:, 1] = reftDF[[1, 2, 3, 4]].agg(" ".join, axis=1)  # dd/mm/yy/time cols are
    reftDF = reftDF.drop(columns=[2, 3, 4])  # read separately; combine into one

    if len(reftDF.columns) != 6:
        log.warning("Raw REFT file for {ssscc} has the incorrect number of columns. Check .CAP file.")
        #   Code will break below

    columns = OrderedDict(  # having this as a dict streamlines next steps
        [
            ("index_memory", int),
            ("datetime", object),
            ("btl_fire_num", int),
            ("diff", int),
            ("raw_value", float),
            ("T90", float),
        ]
    )
    reftDF.columns = list(columns.keys())  # name columns
    reftDF = reftDF.astype(columns)  # force dtypes
    reftDF["datetime"] = pd.to_datetime(reftDF['datetime'], format='%d %b %Y %H:%M:%S') #   Datetime checks, otherwise col is not used

    #   Check contents of the file for stuff that the analyst should double-check (don't immediately flag)
    if any(((reftDF['datetime'] < pd.Timestamp('1985-01-01')) | (reftDF['datetime'] > pd.Timestamp.today()))):
        #   Check timestamps
        log.warning(f"Raw REFT file for {ssscc} contains erroneous timestamps. Check instrument configuration.")
    if any((reftDF['raw_value'] < 50000.0) | (reftDF['raw_value'] > 1000000)):
        #   Check frequencies
        log.warning(f"Raw REFT file reports potentially erroneous frequencies for SSSCC {ssscc}.")
    if any((reftDF['T90'] < -2) | (reftDF['T90'] > 40)):
        #   Check calibration coeffs based on frequencies
        log.warning(f"Raw REFT file reports potentially erroneous temperatures for SSSCC {ssscc}. Check instrument frequencies and calibration coefficients.")
    if any((reftDF['btl_fire_num'] == 1) & (reftDF['btl_fire_num'].shift(1) > 1)):
        #   Check bottle column incrementations
        log.warning(f"Raw REFT file has bottle reset to 1 during {ssscc}. Check .CAP file for tests or other casts.")
    elif reftDF["btl_fire_num"].duplicated().any():
        #   Give warnings if there are any duplicates (file not broken up or deck test bottles not removed)
        log.warning(f"Raw REFT for {ssscc} contains duplicate bottle numbers. Check .CAP file.")

    # assign initial flags (large "diff" = unstable reading, flag questionable)
    reftDF["REFTMP_FLAG_W"] = 2
    reftDF.loc[reftDF["diff"].abs() >= 3000, "REFTMP_FLAG_W"] = 3

    if any(reftDF["REFTMP_FLAG_W"] > 2):    #   Tell the user iteratively in the logs
        for idx, row in reftDF.loc[reftDF["REFTMP_FLAG_W"] == 3].iterrows():
            bad_point = row["index_memory"]
            log.info(f"Measurement {bad_point} flagged questionable in SSSCC {ssscc}")

    # add in STNNBR, CASTNO columns
    # string prob better for other sta/cast formats (names, letters, etc.)
    if len(ssscc) > 5:
        log.warning(f"Length of {ssscc} name exceeds 5. Assigning STNNBR to {ssscc[0:3]}, CASTNO to {ssscc[3:5]}")
    reftDF["STNNBR"] = ssscc[0:3]
    reftDF["CASTNO"] = ssscc[3:5]
    return reftDF


def proc_reft(casts, raw_dir, parsed_dir, cnv_dir):
    """
    This will parse and process all reft data every time it runs, in order to
    make sure any changes or corrections to raw data are picked up and the
    parsed and processed files remain in sync.

    Parameters
    ----------
    casts
    raw_dir
    parsed_dir
    cnv_dir

    Returns
    -------

    """
    parse_sbe35(casts, raw_dir, parsed_dir)

    outdir = validate_dir(cnv_dir, create=True)
    for cast_id in casts:
        infile = Path(parsed_dir, '%s.csv' % cast_id)
        fname = Path(outdir, '%s_reft.csv' % cast_id)

        if infile.exists():
            reft_df = reft_loader(infile)
            reft_df.to_csv(fname, index=False)
        else:
            log.warning("reft file for cast %s cannot be found. Skipping..." % cast_id)



def process_reft(ssscc_list, reft_dir=cfg.dirs["reft"]):
    """
    SBE35 reference thermometer processing function. Load in .cap files for given
    station/cast list, perform basic flagging, and export to .csv files.

    Parameters
    -------
    ssscc_list : list of str
        List of stations to process
    reft_dir : str, optional
        Path to folder containing raw salt files (defaults to data/reft/)

    """
    for ssscc in ssscc_list:
        if Path(reft_dir + ssscc + "_reft.csv").exists():
            try:
                reftDF = _reft_loader(ssscc, reft_dir)
                reftDF.to_csv(reft_dir + ssscc + "_reft.csv", index=False)
            except FileNotFoundError:
                log.warning(
                    "refT file for cast " + ssscc + " does not exist... skipping"
                )
                continue


