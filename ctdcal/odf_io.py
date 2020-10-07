import csv
import datetime as dt
import sys
from collections import OrderedDict
from pathlib import Path

import config as cfg
import gsw
import numpy as np
import pandas as pd


def salts_time_indexer(time_string):
    '''Take in string time, then change it to seconds relative to a 0 point.'''
    try:
        output = dt.datetime.strptime(time_string, f'%H:%M:%S')
    except ValueError:
        output = dt.datetime.strptime(time_string, f'%H:%M:%S*')
    output = output - dt.datetime(1900, 1, 1, 0, 0, 0)
    return output.seconds

def load_short_salts_2017(filepath):
    '''Load a salts file in truncated form and prep for pre-calibration. c. Sept 2017
    Follows normal form.

    Input:
    filepath - a string or path object to the file
    Output:
    dataframe - df following normal form
    '''

    names = [
    'STNNBR',
    'CASTNO',
    'autosal_box_num',
    'bath_temp',
    'cr_average',
    'SAMPNO',
    'standby_number',
    'end_time'
    ]
    usecols = [0,1,2,3,4,5,6,8]

    dtype = {
        'STNNBR':np.int32,
        'CASTNO':np.int32,
        'autosal_box_num':np.int32,
        'bath_temp':np.int32,
        'cr_average':np.float64,
        'SAMPNO':object,
        'standby_number':np.int32
    }
    converters = {
        'end_time': salts_time_indexer
    }

    df = pd.read_csv(filepath, usecols = usecols, names = names, skiprows = [0],
                     delim_whitespace = True, dtype = dtype, converters = converters)

    return df

def salts_create_index_time(df):
    '''Deal with runs that go through midnight, but end less than 24 hours later.'''
    df['index_time'] = df['end_time'] - df['end_time'].iloc[0]
    #make a boolean array of 0s and 1s, then multiply the array by the number to get the number you want, then add together
    df["index_time"] += (df["index_time"] < 0) * (3600*24)
    return df

def create_salt_std_time_array(df):
    '''Create a preformatted array to be passed to np.polyfit for x, or time'''
    df = df.groupby(['SAMPNO'])
    df = df.get_group('worm')
    return df['index_time']

def create_salt_std_cr_array(df):
    '''Create a preformatted array to be passed to np.polyfit for y, or cr'''
    df = df.groupby(['SAMPNO'])
    df = df.get_group('worm')
    return df['cr_average']

def autosal_drift_calibration(df):
    '''Calculate a linear polynomial fit for autosal drift between start and end of runs.'''
    p = np.polyfit(create_salt_std_time_array(df), create_salt_std_cr_array(df), 1)
    p[1] = 0
    p = np.poly1d(p)
    return p

def autosal_drift_fit(df):
    '''Take in the base io dataframe, then compute calibration, then apply fit to data'''
    df = salts_create_index_time(df)
    p = autosal_drift_calibration(df)
    df['cr_average_fitted'] = p(df['index_time']) + df['cr_average']
    df = df[df['SAMPNO'] != 'worm']
    return df

def compute_salinity(df):
    df['SALNTY'] = gsw.SP_salinometer(df['cr_average_fitted']/2, df['bath_temp'])
    return df

def formatted_salt_file(df):
    df = df[['STNNBR', 'CASTNO', 'SAMPNO', 'SALNTY']]
    return df


def _salt_loader(ssscc, salt_dir):
    """
    Load raw file into salt and reference DataFrames.
    """
    saltpath = salt_dir + ssscc  # salt files have no extension
    with open(saltpath, newline="") as f:
        saltF = csv.reader(
            f, delimiter=" ", quoting=csv.QUOTE_NONE, skipinitialspace="True"
        )
        saltArray = []
        for row in saltF:
            saltArray.append(row)
        del saltArray[0]  # remove header

    header = OrderedDict(  # having this as a dict streamlines next steps
        [
            ("STNNBR", int),
            ("CASTNO", int),
            ("SAMPNO", int),
            ("BathTEMP", int),
            ("CRavg", float),
            ("autosalSAMPNO", int),
            ("Unknown", int),
            ("StartTime", object),
            ("EndTime", object),
            ("Attempts", int),
        ]
    )
    saltDF = pd.DataFrame.from_records(saltArray)

    # add as many "Reading#"s as needed
    for ii in range(0, len(saltDF.columns) - len(header)):
        header["Reading{}".format(ii + 1)] = float
    saltDF.columns = list(header.keys())  # name columns

    # TODO: check autosalSAMPNO against SAMPNO for mismatches?
    # TODO: handling for re-samples?

    # add time (in seconds) needed for autosal drift removal step
    saltDF["IndexTime"] = pd.to_datetime(saltDF["EndTime"])
    saltDF["IndexTime"] = (saltDF["IndexTime"] - saltDF["IndexTime"].iloc[0]).dt.seconds
    saltDF["IndexTime"] += (saltDF["IndexTime"] < 0) * (3600 * 24)  # fix overnight runs

    refDF = saltDF.loc[
        saltDF["autosalSAMPNO"] == "worm", ["IndexTime", "CRavg"]
    ].astype(float)
    saltDF = saltDF[saltDF["autosalSAMPNO"] != "worm"].astype(header)  # force dtypes

    return saltDF, refDF


def remove_autosal_drift(saltDF, refDF):
    """Calculate linear CR drift between reference values"""
    diff = refDF.diff(axis="index").dropna()
    time_coef = (diff["CRavg"] / diff["IndexTime"]).iloc[0]

    saltDF["CRavg"] += saltDF["IndexTime"] * time_coef
    saltDF["CRavg"] = saltDF["CRavg"].round(5)  # match initial precision
    saltDF = saltDF.drop(labels="IndexTime", axis="columns")

    return saltDF


def _salt_exporter(saltDF, outdir=cfg.directory["salt"], stn_col="STNNBR", cast_col="CASTNO"):
    """
    Export salt DataFrame to .csv file. Extra logic is included in the event that
    multiple stations and/or casts are included in a single raw salt file.
    """
    stations = saltDF[stn_col].unique()
    for station in stations:
        stn_salts = saltDF[saltDF[stn_col] == station]
        casts = stn_salts[cast_col].unique()
        for cast in casts:
            stn_cast_salts = stn_salts[stn_salts[cast_col] == cast]
            outfile = (  # format to SSSCC_salts.csv
                outdir + "{0:03}".format(station) + "{0:02}".format(cast) + "_salts.csv"
            )
            if Path(outfile).exists():
                print(outfile + " already exists...skipping")
                continue
            stn_cast_salts.to_csv(outfile, index=False)


def process_salts(ssscc_list, salt_dir=cfg.directory["salt"]):
    # TODO: import salt_dir from a config file
    # TODO: update and include linear drift correction from above
    """
    Master salt processing function. Load in salt files for given station/cast list,
    calculate salinity, and export to .csv files.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to process
    salt_dir : str, optional
        Path to folder containing raw salt files (defaults to data/salt/)

    """
    for ssscc in ssscc_list:
        if not Path(salt_dir + ssscc + "_salts.csv").exists():
            try:
                saltDF, refDF = _salt_loader(ssscc, salt_dir)
            except FileNotFoundError:
                print("Salt file for cast " + ssscc + " does not exist... skipping")
                continue
            saltDF = remove_autosal_drift(saltDF, refDF)
            saltDF["SALNTY"] = gsw.SP_salinometer(
                (saltDF["CRavg"] / 2.0), saltDF["BathTEMP"]
            ).round(4)
            _salt_exporter(saltDF, salt_dir)

    return True


def main(argv):
    '''Example script on how to run functions. Not intended for use.'''
    #argv should be a filename
    df = load_short_salts_2017(f'{argv}')
    df = salts_create_index_time(df)
    p = autosal_drift_calibration(df)
    df = autosal_drift_fit(df)
    df = compute_salinity(df)
    df = formatted_salt_file(df)
    return None

if __name__ == '__main__':
    main(sys.argv[1:])
