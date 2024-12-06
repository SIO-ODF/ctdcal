"""Library to create SBE .btl equivalent files.

Joseph Gum SIO/ODF
Nov 7, 2016
"""

import csv
import logging
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from . import flagging as flagging
from . import get_ctdcal_config
from . import oxy_fitting as oxy_fitting
from ctdcal.common import validate_dir

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

BOTTLE_FIRE_COL = "btl_fire"
BOTTLE_FIRE_NUM_COL = "btl_fire_num"


def _get_bottle_order_from_bl_file(bl_file):
    """
    Helper function to get bottle firing order from the Sea Bird .bl file, which
    is the best source when bottles have been fired non-sequentially.

    :param ssscc: string - station/cast identifier
    :return: pandas dataframe object
    """
    btl_fire_order = [0]
    btl_fire_num = [0]
    with open(bl_file, 'r') as bl:
        for line in bl:
            line = line.split(',')
            try:
                int(line[0])
            except ValueError:
                continue
            btl_fire_order.append(int(line[0]))  # index of first bottle fired
            btl_fire_num.append(int(line[1]))  # number of first bottle fired
    return pd.DataFrame(btl_fire_num, index=btl_fire_order, columns=['btl_fire_num'])

def retrieveBottleDataFromFile(converted_file):
    """
    Retrieve the bottle data from a converted file.
    """
    converted_df = pd.read_pickle(converted_file)

    return retrieveBottleData(converted_df)


def retrieveBottleData(converted_df, bl_file):
    """
    Retrieve the bottle data from a dataframe created from a converted file.

    Looks for changes in the BOTTLE_FIRE_COL column, ready to be averaged in making the CTD bottle file.
    """
    bl = _get_bottle_order_from_bl_file(bl_file)
    if BOTTLE_FIRE_COL in converted_df.columns:
        converted_df[BOTTLE_FIRE_NUM_COL] = (
            (
                (converted_df[BOTTLE_FIRE_COL])
                & (
                    converted_df[BOTTLE_FIRE_COL]
                    != converted_df[BOTTLE_FIRE_COL].shift(1)
                )
            )
            .astype(int)
            .cumsum()
        )
        # converted_df['bottle_fire_num'] = ((converted_df[BOTTLE_FIRE_COL] == False)).astype(int).cumsum()
        converted_df[BOTTLE_FIRE_NUM_COL] = bl.loc[converted_df[BOTTLE_FIRE_NUM_COL]].values
        return converted_df.loc[converted_df[BOTTLE_FIRE_COL]]
        # return converted_df
    else:
        log.error("Bottle fire column:", BOTTLE_FIRE_COL, "not found")
        return pd.DataFrame()  # empty dataframe


def bottle_mean(btl_df):
    """Compute the mean for each bottle from a dataframe."""
    btl_max = int(btl_df[BOTTLE_FIRE_NUM_COL].tail(n=1))
    i = 1
    output = pd.DataFrame()
    while i <= btl_max:
        output = pd.concat(
            (
                output,
                btl_df[btl_df[BOTTLE_FIRE_NUM_COL] == i]
                .mean()
                .to_frame(name=i)
                .transpose(),
            )
        )
        i += 1
    return output


def bottle_median(btl_df):
    """Compute the median for each bottle from a dataframe."""
    btl_max = int(btl_df[BOTTLE_FIRE_NUM_COL].tail(n=1))
    i = 1
    output = pd.DataFrame()
    while i <= btl_max:
        output = pd.concat(
            (
                output,
                btl_df[btl_df[BOTTLE_FIRE_NUM_COL] == i]
                .median()
                .to_frame(name=i)
                .transpose(),
            )
        )
        i += 1
    return output


def _load_btl_data(btl_file, cast_id, cols=None):
    """
    Loads "bottle mean" CTD data from .pkl file. Function will return all data unless
    cols is specified (as a list of column names)
    """

    btl_data = pd.read_pickle(btl_file)
    if cols is not None:
        btl_data = btl_data[cols]
    btl_data["SSSCC"] = cast_id

    return btl_data


def _load_reft_data(reft_file, cast_id, index_name="btl_fire_num"):
    """
    Loads reft_file to dataframe and reindexes to match bottle data dataframe
    """
    reft_data = pd.read_csv(reft_file, usecols=["btl_fire_num", "T90", "REFTMP_FLAG_W"])
    reft_data.set_index(index_name)
    reft_data["SSSCC_TEMP"] = cast_id
    reft_data["REFTMP"] = reft_data["T90"]

    return reft_data


def _load_salt_data(salt_file, cast_id, index_name="SAMPNO"):
    """
    Loads salt_file to dataframe and reindexes to match bottle data dataframe
    """
    salt_data = pd.read_csv(
        salt_file, usecols=["SAMPNO", "SALNTY"]
    )
    salt_data.set_index(index_name)
    salt_data["SSSCC_SALT"] = cast_id
    salt_data.rename(columns={"SAMPNO": "SAMPNO_SALT"}, inplace=True)

    return salt_data


def _add_btl_bottom_data(df, cast, datadir, lat_col="LATITUDE", lon_col="LONGITUDE", decimals=4):
    """
    Adds lat/lon, date, and time to dataframe based on the values in the bottom_bottle_details.csv
    """
    bbdfile = Path(datadir, 'logs', 'bottom_bottle_details.csv')
    cast_details = pd.read_csv(bbdfile, dtype={'SSSCC': str})
    cast_details = cast_details[cast_details["SSSCC"] == cast]
    # df[lat_col] = np.round(cast_details["latitude"].iat[0], decimals)
    # df[lon_col] = np.round(cast_details["longitude"].iat[0], decimals)
    df[lat_col] = cast_details["latitude"].iat[0]
    df[lon_col] = cast_details["longitude"].iat[0]

    ts = pd.to_datetime(cast_details["bottom_time"].iat[0], unit="s")
    date = ts.strftime("%Y%m%d")
    hour = ts.strftime("%H%M")
    df["DATE"] = date
    df["TIME"] = hour
    return df


def load_all_btl_files(ssscc_list, datadir, inst, reft, salt, oxy, cols=None):
    """
    Load bottle and secondary (e.g. reference temperature, bottle salts, bottle oxygen)
    files for station/cast list and merge into a dataframe.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to load
    cols : list of str, optional
        Subset of columns to load, defaults to loading all

    Returns
    -------
    df_data_all : DataFrame
        Merged dataframe containing all loaded data

    """
    df_data_all = pd.DataFrame()

    for ssscc in ssscc_list:
        log.info("Loading BTL data for station: " + ssscc + "...")
        btl_file = Path(datadir, 'bottle', inst, "%s_btl_mean.pkl" % ssscc)
        btl_data = _load_btl_data(btl_file, ssscc, cols)

        ### load REFT data
        reft_file = Path(datadir, 'converted', reft, "%s_reft.csv" % ssscc)
        try:
            reft_data = _load_reft_data(reft_file, ssscc)
            if len(reft_data) > 36:
                log.error(f"len(reft_data) > 36 for {ssscc}, check reftmp file")
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFT Data Station: "
                + ssscc
                + "...filling with NaNs"
            )
            reft_data = pd.DataFrame(index=btl_data.index, columns=["T90"], dtype=float)
            reft_data["btl_fire_num"] = btl_data["btl_fire_num"].astype(int)
            reft_data["SSSCC_TEMP"] = ssscc

        ### load REFC data
        refc_file = Path(datadir, 'converted', salt, "%s_salt.csv" % ssscc)
        try:
            refc_data = _load_salt_data(refc_file, ssscc, index_name="SAMPNO")
            if len(refc_data) > 36:
                log.error(f"len(refc_data) > 36 for {ssscc}, check autosal file")
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFC Data Station: "
                + ssscc
                + "...filling with NaNs"
            )
            refc_data = pd.DataFrame(
                index=btl_data.index,
                columns=["SALNTY", "BTLCOND"],
                dtype=float,
            )
            refc_data["SAMPNO_SALT"] = btl_data["btl_fire_num"].astype(int)

        ### load OXY data
        oxy_file = Path(datadir, 'converted', oxy, '%s_oxygen.csv' % ssscc)
        if oxy_file.exists():
            # Added to use non-odf oxy data from an excel sheet and a generic parser
            oxy_data = pd.read_csv(oxy_file, usecols=["SAMPNO", "OXYGEN"])
            oxy_data['SSSCC_OXY'] = ssscc
            oxy_data.rename(columns={'SAMPNO': 'BOTTLENO_OXY'}, inplace=True)
        else:
            oxy_file = Path(cfg.dirs["oxygen"] + ssscc)
            try:
                oxy_data, params = oxy_fitting.load_winkler_oxy(oxy_file)
                if len(oxy_data) > 36:
                    log.error(f"len(oxy_data) > 36 for {ssscc}, check oxygen file")
            except FileNotFoundError:
                log.warning(
                    "Missing (or misnamed) REFO Data Station: "
                    + ssscc
                    + "...filling with NaNs"
                )
                oxy_data = pd.DataFrame(
                    index=btl_data.index,
                    columns=[
                        "FLASKNO",
                        "TITR_VOL",
                        "TITR_TEMP",
                        "DRAW_TEMP",
                        "TITR_TIME",
                        "END_VOLTS",
                    ],
                    dtype=float,
                )
                oxy_data["STNNO_OXY"] = ssscc[:3]
                oxy_data["CASTNO_OXY"] = ssscc[3:]
                oxy_data["BOTTLENO_OXY"] = btl_data["btl_fire_num"].astype(int)

        ### clean up dataframe
        # Horizontally concat DFs to have all data in one DF
        # btl_data = pd.merge(btl_data, reft_data, on="btl_fire_num", how="outer")
        btl_data = pd.merge(btl_data, reft_data, on="btl_fire_num", how="left")
        btl_data = pd.merge(
            btl_data,
            refc_data,
            left_on="btl_fire_num",
            right_on="SAMPNO_SALT",
            how="outer",
        )
        btl_data = pd.merge(
            btl_data,
            oxy_data,
            left_on="btl_fire_num",
            right_on="BOTTLENO_OXY",
            how="outer",
        )

        if len(btl_data) > 36:
            log.error(f"len(btl_data) for {ssscc} > 36, check bottle file")

        # Add bottom of cast information (date,time,lat,lon,etc.)
        btl_data = _add_btl_bottom_data(btl_data, ssscc, datadir)

        # Merge cast into df_data_all
        try:
            df_data_all = pd.concat([df_data_all, btl_data], sort=False)
        except AssertionError:
            raise AssertionError(
                "Columns of " + ssscc + " do not match those of previous columns"
            )
        # print("* Finished BTL data station: " + ssscc + " *")

    # Drop duplicated columns generated by concatenation
    df_data_all = df_data_all.loc[:, ~df_data_all.columns.duplicated()]

    df_data_all["master_index"] = range(len(df_data_all))

    return df_data_all


def _reft_loader(ssscc, search_dir):
    """
    Loads SBE35.cap files and assembles into a dataframe.

    Parameters
    ----------
    ssscc : str
        Station to load .CAP file for
    search_dir : str, pathlib Path
        Path to the reft folder 

    Returns
    -------
    reftDF : DataFrame
        DataFrame of .CAP file with headers

    """
    # semi-flexible search for reft file (in the form of *ssscc.cap OR .txt)
    try:
        reft_path = sorted(Path(search_dir).glob(f"*{ssscc}.[ct][ax][pt]"))[0]
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


def process_reft(ssscc_list, data_dir, inst):
    """
    SBE35 reference thermometer processing function. Load in .cap files for given
    station/cast list, perform basic flagging, and export to .csv files.

    Parameters
    -------
    ssscc_list : list of str
        List of stations to process
    raw_dir : str, optional
        Path to folder containing raw salt files (defaults to data/reft/)

    """
    raw_dir = Path(data_dir, 'raw', inst)
    out_dir = validate_dir(Path(data_dir, 'converted', inst), create=True)
    for ssscc in ssscc_list:
        fname = Path(out_dir, "%s_reft.csv" % ssscc)
        if not fname.exists():
            try:
                reftDF = _reft_loader(ssscc, raw_dir)
                reftDF.to_csv(fname, index=False)
            except FileNotFoundError:
                log.warning(
                    "refT file for cast " + ssscc + " does not exist... skipping"
                )
                continue


def add_btlnbr_cols(df, btl_num_col):
    """
    Initialize bottle number column and initialize WOCE bottle flags.

    Parameters
    ----------
    df : Pandas DataFrame
        Bottle DataFrame containing a defined rosette bottle number
    btl_num_col : String
        String of bottle column to be reassigned

    Returns
    -------
    df : Pandas DataFrame
        Bottle DataFrame with BTLNBR and flag columns as type int
    """
    df["BTLNBR"] = df[btl_num_col].astype(int)
    # default to everything being good
    df["BTLNBR_FLAG_W"] = 2
    return df


def load_hy_file(path_to_hyfile):
    """
    Read in an exchange-formatted bottle file as a Pandas DataFrame.

    Inputs
    ------
    path_to_hyfile : String or Path object
        The path to the bottle file.

    Returns
    -------
    df : Pandas DataFrame
        The bottle file without the lead/end rows, comments, or units
    """

    df = pd.read_csv(path_to_hyfile, comment="#", skiprows=[0])
    df = df.drop(df.index[0])  #   Drop the units
    df = df[df["EXPOCODE"] != "END_DATA"]  #   Drop the final row
    return df


def merge_hy1(
    df1,
    df2,
):
    """
    Merges two hy1 files, returning the combined Pandas DataFrame.
    If the hy1 file has not been loaded yet, use load_hy_file.

    Inputs
    -------
    df1 : Pandas DataFrame
        First hy1 file for concatination
    df2 : Pandas DataFrame
        Second hy1 file for concatination

    Returns
    df: Pandas DataFrame
        Merged bottle file as a DataFrame
    """

    if set(df1.columns) != set(df2.columns):
        print("Bottle file columns do not match. Concatenating with NaNs.")

    df = pd.concat([df1, df2], axis=0, ignore_index=True)  #   Staple df2 onto there

    sorting_cols = {"STNNBR", "CASTNO", "SAMPNO"}
    if sorting_cols.issubset(df):
        if df[list(sorting_cols)].isna().any().any():
            print("NaNs found in station/cast/sample number. Check source files.")

        else:
            df = df.sort_values(
                by=["STNNBR", "CASTNO", "SAMPNO"],
                ascending=[True, True, False],
                ignore_index=True,
            )
    return df


def export_report_data(df, datadir, inst):
    """
    Write out the data used for report generation as a csv.

    Params
    ------
    df : Pandas DataFrame
        Fit bottle data

    """
    outdir = validate_dir(Path(datadir, 'report', inst), create=True)
    df["CAST"] = df["SSSCC"]
    df["CTDPRS"] = df["CTDPRS"].round(1)
    cruise_report_cols = [
        "CAST",
        "CTDPRS",
        "CTDTMP1",
        "CTDTMP1_FLAG_W",
        "CTDTMP2",
        "CTDTMP2_FLAG_W",
        "REFTMP",
        "CTDCOND1",
        "CTDCOND1_FLAG_W",
        "CTDCOND2",
        "CTDCOND2_FLAG_W",
        "BTLCOND",
        "CTDSAL",
        "CTDSAL_FLAG_W",
        "SALNTY",
        "CTDOXY",
        "CTDOXY_FLAG_W",
        # "CTDRINKO",
        # "CTDRINKO_FLAG_W",
        "OXYGEN",
    ]

    # add in missing flags
    df["CTDTMP1_FLAG_W"] = flagging.by_residual(
        df["CTDTMP1"], df["REFTMP"], df["CTDPRS"]
    )
    df["CTDTMP2_FLAG_W"] = flagging.by_residual(
        df["CTDTMP1"], df["REFTMP"], df["CTDPRS"]
    )
    df["CTDCOND1_FLAG_W"] = flagging.by_residual(
        df["CTDCOND1"], df["BTLCOND"], df["CTDPRS"]
    )
    df["CTDCOND2_FLAG_W"] = flagging.by_residual(
        df["CTDCOND2"], df["BTLCOND"], df["CTDPRS"]
    )
    df["CTDOXY_FLAG_W"] = flagging.by_percent_diff(df["CTDOXY"], df["OXYGEN"])
    # df["CTDRINKO_FLAG_W"] = flagging.by_percent_diff(df["CTDRINKO"], df["OXYGEN"])

    df[cruise_report_cols].to_csv(Path(outdir, "report_data.csv"), index=False)

    return


def export_hy1(df, datadir, inst, org="ODF"):
    """
    Write out the exchange-lite formatted hy1 bottle file.

    Params
    ------
    df : Pandas DataFrame
        Fit bottle data
    out_dir = String or Path object, optional
        The path for where to write the hy1 file
    org : String, optional
        The organization or group used to determine subroutines

    """
    log.info("Exporting bottle file")
    btl_data = df.copy()
    now = datetime.now()
    file_datetime = now.strftime("%Y%m%d")

    btl_columns = {
        "EXPOCODE": "",
        "SECT_ID": "",
        # "STNNBR": "",
        "CASTNO": "",
        "SAMPNO": "",
        "BTLNBR": "",
        "BTLNBR_FLAG_W": "",
        "DATE": "",
        "TIME": "",
        "LATITUDE": "",
        "LONGITUDE": "",
        "DEPTH": "METERS",
        "CTDPRS": "DBAR",
        "CTDTMP": "ITS-90",
        "CTDSAL": "PSS-78",
        "CTDSAL_FLAG_W": "",
        "SALNTY": "PSS-78",
        "SALNTY_FLAG_W": "",
        # "CTDOXY": "UMOL/KG",
        # "CTDOXY_FLAG_W": "",
        # "CTDRINKO": "UMOL/KG",
        # "CTDRINKO_FLAG_W": "",
        "CTDOXY": "UMOL/KG",
        "CTDOXY_FLAG_W": "",
        "OXYGEN": "UMOL/KG",
        "OXYGEN_FLAG_W": "",
        "REFTMP": "ITS-90",
        "REFTMP_FLAG_W": "",
    }

    # rename outputs as defined in user_settings.yaml
    for param, attrs in cfg.ctd_outputs.items():
        if param not in btl_data.columns:
            btl_data.rename(columns={attrs["sensor"]: param}, inplace=True)

    btl_data["EXPOCODE"] = cfg.expocode
    btl_data["SECT_ID"] = cfg.section_id
    # btl_data["STNNBR"] = [int(x[0:3]) for x in btl_data["SSSCC"]]
    btl_data["CASTNO"] = btl_data["SSSCC"].astype(str)
    btl_data["SAMPNO"] = btl_data["btl_fire_num"].astype(int)
    btl_data = add_btlnbr_cols(btl_data, btl_num_col="btl_fire_num")

    # sort by decreasing sample number (increasing pressure) and reindex
    btl_data = btl_data.sort_values(
        by=["CASTNO", "SAMPNO"], ascending=[True, True], ignore_index=True
    )

    # switch oxygen primary sensor to rinko
    # btl_data["CTDOXY"] = btl_data.loc[:, "CTDRINKO"]
    # btl_data["CTDOXY_FLAG_W"] = btl_data.loc[:, "CTDRINKO_FLAG_W"]

    # round data
    # for col in ["CTDTMP", "CTDSAL", "SALNTY", "REFTMP"]:
    #     btl_data[col] = btl_data[col].round(4)
    # for col in ["CTDPRS", "CTDOXY", "OXYGEN"]:
    #     btl_data[col] = btl_data[col].round(1)

    # add depth
    depth_df = pd.read_csv(
        Path(datadir, "logs", "depth_log.csv"), dtype={"SSSCC": str}, na_values=-999
    ).dropna()
    manual_depth_df = pd.read_csv(
        Path(datadir, "logs", "manual_depth_log.csv"), dtype={"SSSCC": str}
    )
    full_depth_df = pd.concat([depth_df, manual_depth_df])
    full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)
    btl_data["DEPTH"] = -999
    for index, row in full_depth_df.iterrows():
        btl_data.loc[btl_data["SSSCC"] == row["SSSCC"], "DEPTH"] = int(row["DEPTH"])

    # deal with nans
    btl_data["REFTMP_FLAG_W"] = flagging.nan_values(
        btl_data["REFTMP_FLAG_W"], old_flags=btl_data["REFTMP_FLAG_W"]
    )
    btl_data = btl_data.where(~btl_data.isnull(), -999)

    # check columns
    try:
        btl_data[btl_columns.keys()]
        # this is lazy, do better
    except KeyError as err:
        log.info("Column names not configured properly... attempting to correct")
        bad_cols = err.args[0].split("'")[1::2]  # every other str is a column name
        for col in bad_cols:
            if col.endswith("FLAG_W"):
                log.warning(col + " missing, flagging with 9s")
                btl_data[col] = 9
            else:
                log.warning(col + " missing, filling with -999s")
                btl_data[col] = -999

    btl_data = btl_data[btl_columns.keys()]
    time_stamp = file_datetime + org
    with open(Path(datadir, "export", inst, "%s_hy1.csv" % cfg.expocode), mode="w+") as f:
        f.write("BOTTLE, %s\n" % (time_stamp))
        f.write(",".join(btl_columns.keys()) + "\n")
        f.write(",".join(btl_columns.values()) + "\n")
        btl_data.to_csv(f, header=False, index=False)
        f.write("\n" + "END_DATA")

    return
