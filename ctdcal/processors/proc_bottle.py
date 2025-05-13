"""
Transforming processed bottle data into derivative forms (binned ctd bottle data,
merged water sample data)
"""
import logging
from pathlib import Path

import pandas as pd

from ctdcal import get_ctdcal_config
from ctdcal.common import validate_dir
from ctdcal.processors.proc_oxy_odf import load_winkler_oxy


cfg = get_ctdcal_config()  # legacy config
log = logging.getLogger(__name__)


def make_btl_files(casts, raw_dir, btl_dir, cnv_dir):
    log.info("Generating btl_mean.pkl files")
    # validate required directories
    raw_dir = validate_dir(Path(raw_dir))
    btl_dir = validate_dir(Path(btl_dir), create=True)
    cnv_dir = validate_dir(Path(cnv_dir))
    # log_dir = validate_dir(Path(log_dir))

    for cast_id in casts:
        btl_file = Path(btl_dir, '%s_btl_mean.pkl' % cast_id)
        if not btl_file.exists():
            cnv_file = Path(cnv_dir, '%s.pkl' % cast_id)
            cnv_df = pd.read_pickle(cnv_file)
            firing_order = get_bottle_order_from_bl_file(cast_id, raw_dir)
            bottle_df = get_bottle_data(cnv_df, firing_order)
            mean_df = bottle_df.groupby('btl_fire_num', as_index=False).mean()
            mean_df['cast_id'] = cast_id
            mean_df.to_pickle(btl_file)
    return True


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
    log.warning("Use of make_btl_mean() is deprecated. Use make_btl_files() instead.")
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


def get_bottle_data(cnv, firing_order):
    if 'btl_fire' in cnv.columns:
        # find all rows associated with bottle fires
        cnv['btl_fire_num'] = (
            (
                (cnv['btl_fire'])
                & (
                    cnv['btl_fire']
                    != cnv['btl_fire'].shift(1)
                )
            )
            .astype(int)
            .cumsum()
        )
        # label rows according to actual firing order
        cnv['btl_fire_num'] = firing_order.loc[cnv['btl_fire_num']].values
        return cnv.loc[cnv['btl_fire']]
    else:
        log.error("Bottle fire column not found")
        return pd.DataFrame()  # empty dataframe


def get_bottle_order_from_bl_file(cast_id, raw_dir):
    """
    Helper function to get bottle firing order from the Sea Bird .bl file, which
    is the best source when bottles have been fired non-sequentially.

    :param cast_id: string - cast name

    :return: pandas dataframe object
    """
    f = Path(raw_dir, '%s.bl' % cast_id)
    btl_fire_order = [0]
    btl_fire_num = [0]
    with open(f, 'r') as bl:
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
    log.warning("Use of retrieveBottleData() is deprecated. Use get_btl_data() instead.")
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


def get_bottom_bottle_data(
        btl_data,
        report_dir,
        cols=None,
        export=False
):
    if cols is None:
        cols = {'cast_id': 'SSSCC', 'date': 'DATE', 'time': 'TIME', 'lat': 'LATITUDE', 'lon': 'LONGITUDE'}
    # datetime_col = "nmea_datetime"
    # ## 2025-04-23 doing that clobbers casts without nmea datetime
    datetime_col = "scan_datetime"
    # TODO: make this determination cast-by-cast because maybe this wasn't enabled the whole time ; )
    if datetime_col not in btl_data.columns:
        log.debug(
                f"'{datetime_col}' not found in DataFrame - using 'scan_datetime'"
        )
        datetime_col = "scan_datetime"

    # filter bottom bottle data
    bottom_bottles = btl_data.groupby('cast_id', as_index=False)[[datetime_col, 'GPSLAT', 'GPSLON']].first()

    # convert timestamps
    bottom_bottles[datetime_col] = pd.to_datetime(bottom_bottles[datetime_col], unit="s")
    bottom_bottles[cols['date']] = bottom_bottles[datetime_col].dt.strftime("%Y%m%d")
    bottom_bottles[cols['time']] = bottom_bottles[datetime_col].dt.strftime("%H%M")

    # rename to user cols
    bottom_bottles.rename(
            columns={'cast_id': cols['cast_id'],
                     'GPSLAT': cols['lat'],
                     'GPSLON': cols['lon'],
                     },
            inplace=True,
    )
    bottom_bottles.drop(columns=datetime_col, inplace=True)

    # export bottom bottle time/lat/lon info
    if export is True:
        bottom_bottle_file = Path(report_dir, 'bottom_bottle_details.csv')
        bottom_bottles.to_csv(bottom_bottle_file, header=True, index=False)

    return bottom_bottles




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


def load_btl_all(casts, btl_dir, reft_dir, salt_dir, oxy_dir):
    btl_dir = Path(btl_dir)
    reft_dir = Path(reft_dir)
    salt_dir = Path(salt_dir)

    reft_data_all = pd.DataFrame()
    refc_data_all = pd.DataFrame()
    oxy_data_all = pd.DataFrame()
    btl_data_all = pd.DataFrame()

    for cast_id in casts:
        log.info("Loading BTL data for station: " + cast_id + "...")
        btl_file = Path(btl_dir, "%s_btl_mean.pkl" % cast_id)
        btl_data = pd.read_pickle(btl_file)

        ### load REFT data
        reft_file = Path(reft_dir, "%s_reft.csv" % cast_id)
        try:
            reft_data = pd.read_csv(reft_file, dtype={'cast_id': str}, na_values='NaN')
            if len(reft_data) > 36:
                log.error(f"len(reft_data) > 36 for {cast_id}, check reftmp file")
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFT Data Station: "
                + cast_id
                + "...filling with NaNs"
            )
            reft_data = pd.DataFrame()
            # reft_data = pd.DataFrame(index=btl_data.index, columns=["REFTMP"], dtype=float)
            # reft_data["btl_fire_num"] = btl_data["btl_fire_num"].astype(int)
            # reft_data["cast_id"] = cast_id

            # concatenate with all reft
        reft_data_all = pd.concat([reft_data_all, reft_data], axis=0, sort=False)

        ### load REFC data
        refc_file = Path(salt_dir, "%s_salts.csv" % cast_id)
        try:
            refc_data = pd.read_csv(refc_file, dtype={'cast_id': str}, na_values='NaN')
            if len(refc_data) > 36:
                log.error(f"len(refc_data) > 36 for {cast_id}, check autosal file")
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFC Data Station: "
                + cast_id
                + "...filling with NaNs"
            )
            refc_data = pd.DataFrame()
            # refc_data = pd.DataFrame(
            #     index=btl_data.index,
            #     columns=["CRavg", "BathTEMP", "BTLCOND"],
            #     dtype=float,
            # )
            # refc_data["btl_fire_num"] = btl_data["btl_fire_num"].astype(int)
            # refc_data['cast_id'] = cast_id

        # concatenate with all refc
        refc_data_all = pd.concat([refc_data_all, refc_data], axis=0, sort=False)

        ### load OXY data
        oxy_file = Path(oxy_dir, '%s_oxy.csv' % cast_id)
        try:
            oxy_data = pd.read_csv(oxy_file, dtype={'cast_id': str, 'FLASKNO': str}, na_values='NaN')
            if len(oxy_data) > 36:
                log.error(f"len(oxy_data) > 36 for {cast_id}, check oxygen file")
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFO Data Station: "
                + cast_id
                + "...filling with NaNs"
            )
            oxy_data = pd.DataFrame()
            # oxy_data = pd.DataFrame(
            #     index=btl_data.index,
            #     columns=[
            #         "FLASKNO",
            #         "TITR_VOL",
            #         "TITR_TEMP",
            #         "DRAW_TEMP",
            #         "TITR_TIME",
            #         "END_VOLTS",
            #     ],
            #     dtype=float,
            # )
            # oxy_data["STNNO_OXY"] = cast_id[:3]
            # oxy_data["CASTNO_OXY"] = cast_id[3:]
            # oxy_data["btl_fire_num"] = btl_data["btl_fire_num"].astype(int)
            # oxy_data['cast_id'] = cast_id

        # concatenate with all oxy
        oxy_data_all = pd.concat([oxy_data_all, oxy_data], axis=0, sort=False)

        # concatenate with all btl
        btl_data_all = pd.concat([btl_data_all, btl_data], axis=0, sort=False)

    ### build up dataframe
    # Horizontally concat DFs to have all data in one DF
    btl_data_all = pd.merge(btl_data_all, reft_data_all, on=['cast_id', 'btl_fire_num'], how='left')
    btl_data_all = pd.merge(btl_data_all, refc_data_all, on=['cast_id', 'btl_fire_num'], how='left')
    btl_data_all = pd.merge(btl_data_all, oxy_data_all, on=['cast_id', 'btl_fire_num'], how='left')

    return btl_data_all


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
    log.warning("Use of load_all_btl_files() is deprecated. Use load_btl_all() instead.")
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
