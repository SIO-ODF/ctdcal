"""
Transforming processed bottle data into derivative forms (binned ctd bottle data,
merged water sample data)
"""
import logging
from pathlib import Path

import pandas as pd

from ctdcal import get_ctdcal_config
from ctdcal.processors.proc_oxy_odf import load_winkler_oxy

# from ctdcal.process_bottle import retrieveBottleData, bottle_mean

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def make_btl_mean(ssscc_list):
    """
    Create "bottle mean" files from continuous CTD data averaged at the bottle stops.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to convert

    Returns
    -------
    boolean
        bottle averaging of mean has finished successfully
    """
    log.info("Generating btl_mean.pkl files")
    for ssscc in ssscc_list:
        if not Path(cfg.dirs["bottle"] + ssscc + "_btl_mean.pkl").exists():
            imported_df = pd.read_pickle(cfg.dirs["converted"] + ssscc + ".pkl")
            bottle_df = retrieveBottleData(imported_df)
            mean_df = bottle_mean(bottle_df)

            # export bottom bottle time/lat/lon info
            fname = cfg.dirs["logs"] + "bottom_bottle_details.csv"
            datetime_col = "nmea_datetime"
            if datetime_col not in mean_df.columns:
                log.debug(
                    f"'{datetime_col}' not found in DataFrame - using 'scan_datetime'"
                )
                datetime_col = "scan_datetime"

            bot_df = mean_df[[datetime_col, "GPSLAT", "GPSLON"]].head(1)
            bot_df.columns = ["bottom_time", "latitude", "longitude"]
            bot_df.insert(0, "SSSCC", ssscc)
            add_header = not Path(fname).exists()  # add header iff file doesn't exist
            with open(fname, "a") as f:
                bot_df.to_csv(f, mode="a", header=add_header, index=False)

            mean_df.to_pickle(cfg.dirs["bottle"] + ssscc + "_btl_mean.pkl")

    return True


def retrieveBottleDataFromFile(converted_file):
    """
    Retrieve the bottle data from a converted file.
    """
    converted_df = pd.read_pickle(converted_file)

    return retrieveBottleData(converted_df)


# Presumably there's a historical reason that these values are defined as constants, but
# this doesn't seem the best place to do it. If they need to be redefined (maybe to process
# older data?) they should be set in the user config file.
# TODO: move these constants to user config, or go with fixed values
BOTTLE_FIRE_COL = "btl_fire"
BOTTLE_FIRE_NUM_COL = "btl_fire_num"
def retrieveBottleData(converted_df):
    """
    Retrieve the bottle data from a dataframe created from a converted file.

    Looks for changes in the BOTTLE_FIRE_COL column, ready to be averaged in making the CTD bottle file.
    """
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


def _load_btl_data(btl_file, cols=None):
    """
    Loads "bottle mean" CTD data from .pkl file. Function will return all data unless
    cols is specified (as a list of column names)
    """

    btl_data = pd.read_pickle(btl_file)
    if cols is not None:
        btl_data = btl_data[cols]
    btl_data["SSSCC"] = Path(btl_file).stem.split("_")[0]

    return btl_data


def _load_reft_data(reft_file, index_name="btl_fire_num"):
    """
    Loads reft_file to dataframe and reindexes to match bottle data dataframe
    """
    reft_data = pd.read_csv(reft_file, usecols=["btl_fire_num", "T90", "REFTMP_FLAG_W"])
    reft_data.set_index(index_name)
    reft_data["SSSCC_TEMP"] = Path(reft_file).stem.split("_")[0]
    reft_data["REFTMP"] = reft_data["T90"]

    return reft_data


def _load_salt_data(salt_file, index_name="SAMPNO"):
    """
    Loads salt_file to dataframe and reindexes to match bottle data dataframe
    """
    salt_data = pd.read_csv(
        salt_file, usecols=["SAMPNO", "SALNTY", "BathTEMP", "CRavg"]
    )
    salt_data.set_index(index_name)
    salt_data["SSSCC_SALT"] = Path(salt_file).stem.split("_")[0]
    salt_data.rename(columns={"SAMPNO": "SAMPNO_SALT"}, inplace=True)

    return salt_data


def _add_btl_bottom_data(df, cast, lat_col="LATITUDE", lon_col="LONGITUDE", decimals=4):
    """
    Adds lat/lon, date, and time to dataframe based on the values in the bottom_bottle_details.csv
    """
    cast_details = pd.read_csv(
        # cfg.dirs["logs"] + "cast_details.csv", dtype={"SSSCC": str}
        cfg.dirs["logs"] + "bottom_bottle_details.csv",
        dtype={"SSSCC": str},
    )
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


def load_all_btl_files(ssscc_list, cols=None):
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
        btl_file = cfg.dirs["bottle"] + ssscc + "_btl_mean.pkl"
        btl_data = _load_btl_data(btl_file, cols)

        ### load REFT data
        reft_file = cfg.dirs["reft"] + ssscc + "_reft.csv"
        try:
            reft_data = _load_reft_data(reft_file)
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
        refc_file = cfg.dirs["salt"] + ssscc + "_salts.csv"
        try:
            refc_data = _load_salt_data(refc_file, index_name="SAMPNO")
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
                columns=["CRavg", "BathTEMP", "BTLCOND"],
                dtype=float,
            )
            refc_data["SAMPNO_SALT"] = btl_data["btl_fire_num"].astype(int)

        ### load OXY data
        oxy_file = Path(cfg.dirs["oxygen"] + ssscc)
        try:
            oxy_data, params = load_winkler_oxy(oxy_file)
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
        btl_data = pd.merge(btl_data, reft_data, on="btl_fire_num", how="outer")
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
        btl_data = _add_btl_bottom_data(btl_data, ssscc)

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
