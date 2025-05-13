"""
Module for legacy SBE conversion and processing.
* Most or all of this module will be deprecated in CTDCAL 1.0

A module for handling SeaBird raw .HEX files, including the generation of bottle-extractions,
downcast isolation, and SBE3/4C handling.
"""
import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from ctdcal.common import validate_dir
from ctdcal.processors import functions_ctd as sbe_eq, sbe_reader as sbe_rd
from ctdcal.processors import functions_oxy as oxy_eq
from ctdcal.processors import functions_aux as aux_eq
from ctdcal import get_ctdcal_config

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

# lookup table for sensor data
# DOUBLE CHECK TYPE IS CORRECT #
short_lookup = {
    "55": {
        "short_name": "CTDTMP",
        "long_name": "SBE 3+ Temperature",
        "units": "ITS-90",
        "type": "float64",
    },
    "45": {
        "short_name": "CTDPRS",
        "long_name": "SBE 9+ Pressure",
        "units": "DBAR",
        "type": "float64",
    },
    "3": {
        "short_name": "CTDCOND",
        "long_name": "SBE 4 Conductivity",
        "units": "MSPCM",
        "type": "float64",
    },
    "38": {
        "short_name": "CTDOXY",
        "long_name": "SBE 43 Oxygen",
        "units": "MLPL",
        "type": "float64",
    },
    # '38':{'short_name': 'CTDOXYVOLTS', 'long_name':'SBE 43 Oxygen Volts', 'units': '0-5VDC', 'type':'float64'},
    "11": {
        "short_name": "FLUOR",
        "long_name": "Seapoint Fluorometer",
        "units": "0-5VDC",
        "type": "float64",
    },
    "27": {"short_name": "FREE", "long_name": "empty", "units": "NA", "type": "NA"},
    "0": {
        "short_name": "ALT",
        "long_name": "Altitude",
        "units": "M",
        "type": "float64",
    },
    "71": {
        "short_name": "CTDXMISS",
        "long_name": "CStar",
        "units": "0-5VDC",
        "type": "float64",
    },
    "61": {
        "short_name": "U_DEF_poly",
        "long_name": "user defined, polynomial",
        "units": "0-5VDC",
        "type": "float64",
    },
    "80": {
        "short_name": "U_DEF_e",
        "long_name": "user defined, exponential",
        "units": "0-5VDC",
        "type": "float64",
    },
    "1000": {
        "short_name": "CTDSAL",
        "long_name": "Salinity (C1 T1)",
        "units": "PSU",
        "type": "float64",
    },
    "20": {
        "short_name": "CTDFLUOR",
        "long_name": "WetlabECO_AFL_FL_Sensor",
        "units": "0-5VDC",
        "type": "float64",
    },
    "42": {
        "short_name": "PAR",
        "long_name": "PAR/Irradiance, Biospherical/Licor",
        "units": "0-5VDC",
        "type": "float64",
    },
    "51": {
        "short_name": "REF_PAR",
        "long_name": "Surface PAR/Irradiance, Biospherical/Licor",
        "units": "0-5VDC",
        "type": "float64",
    },
    "70": {
        "short_name": "CTDBACKSCATTER",
        "long_name": "WetlabECO_BB_Sensor",
        "units": "0-5VDC",
        "type": "float64",
    },
    "13": {
        "short_name": "CTDFLUOR",
        "long_name": "FluoroSeatechWetlabsFLF_Sensor",
        "units": "0-5VDC",
        "type": "float64",
    },
    "67": {
        "short_name": "TURBIDITY",
        "long_name": "Turbidity Meter",
        "units": "0-5VDC",
        "type": "float64",
    },
    "19": {
        "short_name": "FLUOR_CDOM",
        "long_name": "FluoroWetlabCDOM_Sensor",
        "units": "0-5VDC",
        "type": "float64",
    },
}


def hex_to_ctd(ssscc_list, rawdir, outdir):
    """
    Convert raw CTD data and export to .pkl files.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to convert

    Returns
    -------

    """
    log.info("Converting .hex files")
    for ssscc in ssscc_list:
        validate_dir(rawdir, create=True)
        validate_dir(outdir, create=True)
        outfile = Path(outdir, '%s.pkl' % ssscc)
        if not outfile.exists():
            hexFile = Path(rawdir, '%s.hex' % ssscc)
            xmlconFile = Path(rawdir, '%s.XMLCON' % ssscc)
            sbeReader = sbe_rd.SBEReader.from_paths(hexFile, xmlconFile)
            converted_df = convertFromSBEReader(sbeReader, ssscc)
            converted_df.to_pickle(outfile)

    return True

def convertFromSBEReader(sbeReader, ssscc):
    """Handler to convert engineering data to sci units automatically.
    Takes SBEReader object that is already connected to the .hex and .XMLCON files.
    """

    # Retrieve parsed scans and convert to dataframe
    rawData = sbeReader.parsed_scans
    raw_df = pd.DataFrame(rawData).apply(pd.to_numeric, errors="ignore")
    raw_df.index.name = "index"

    # Metadata needs to be processed seperately and then joined with the converted data
    log.info(f"Building metadata dataframe for {ssscc}")
    metaArray = [line.split(",") for line in sbeReader._parse_scans_meta().tolist()]
    meta_cols, meta_dtypes = sbeReader._breakdown_header()
    meta_df = pd.DataFrame(metaArray)
    meta_df.columns = meta_cols
    meta_df.index.name = "index"

    for col, dtype in zip(meta_cols, meta_dtypes):
        if dtype != "bool_":
            meta_df[col] = meta_df[col].astype(dtype)
        else:
            # map from string "pseudo-boolean" to actual boolean values
            meta_df[col] = meta_df[col].map({"True": True, "False": False})

    log.info("Success!")

    t_probe = meta_df["pressure_temp_int"].tolist()  # raw int from Digitquartz T probe

    # Temporary arrays to hold scientific values needed to compute cond/oxy
    t_array, p_array, c_array = [], [], []

    # needs to search sensor dictionary, and compute in order:
    # temp, pressure, cond, salinity, oxygen, all aux.
    # run one loop that builds a queue to determine order of processing, must track which column to pull
    # process queue, store results in seperate arrays for reuse later
    # once queue is empty, attach results together according to format order or xmlcon order - structure to keep track
    queue_metadata = []
    temp_counter = 0
    cond_counter = 0
    oxygen_counter = 0
    u_def_p_counter = 0
    u_def_e_counter = 0
    empty_counter = 0

    # The following are definitions for every key in the dict below:
    #
    # sensor_id = number assigned by SBE for identification in XML
    # list_id = place in XML array by SBE for determining which sensor is which, alternatively channel number (freq+volt)
    # channel_pos = is it the first, second, third, etc sensor of its type in the data file, aux sensors default to 0
    # ranking = data processing ranking - temp first, then pressure, then conductivity, then oxygen, then aux
    # column = column in the raw_df containing the engineering units to be converted to sci units
    # sensor_info = xml sensor info to convert from eng units to sci units

    rawConfig = sbeReader.parsed_config()  # Retrieve Config data

    for list_id, sensor_info in rawConfig["Sensors"].items():
        sensor_id = sensor_info["SensorID"]

        if sensor_id == "55":  # temperature block
            temp_counter += 1
            channel_pos = temp_counter
            ranking = 1

        elif sensor_id == "45":  # pressure block
            channel_pos = ""
            ranking = 2

        elif sensor_id == "3":  # conductivity block
            cond_counter += 1
            channel_pos = cond_counter
            ranking = 3

        elif sensor_id == "38":  # oxygen block
            oxygen_counter += 1
            channel_pos = oxygen_counter
            ranking = 5

        elif sensor_id == "27":  # empty block
            empty_counter += 1
            channel_pos = empty_counter
            ranking = 6

        elif sensor_id == "61":  # user defined (polynomial) block
            u_def_p_counter += 1
            channel_pos = u_def_p_counter
            ranking = 6

        elif sensor_id == "80":  # user defined (exponential) block
            u_def_e_counter += 1
            channel_pos = u_def_e_counter
            ranking = 6

        else:  # auxiliary block
            channel_pos = ""
            ranking = 7

        queue_metadata.append(
            {
                "sensor_id": sensor_id,
                "list_id": list_id,
                "channel_pos": channel_pos,
                "ranking": ranking,
                "column": list_id,
                "sensor_info": sensor_info,
            }
        )

    # Temporary block in order to append basic salinity (calculated from T1/C1) to file
    # If additional salinity is needed (e.g. T2/C2), it'll need a full reworking
    queue_metadata.append(
        {
            "sensor_id": "1000",
            "list_id": 1000,
            "channel_pos": "",
            "ranking": 4,
            "column": "",
            "sensor_info": "",
        }
    )

    # Queue sorting forces it to be in order, so we don't worry about order here
    # Assumes first channel for each sensor is primary for computing following data
    queue_metadata = sorted(queue_metadata, key=lambda sensor: sensor["ranking"])

    # Initialize converted dataframe
    converted_df = pd.DataFrame()

    for meta in queue_metadata:

        col = f"{short_lookup[meta['sensor_id']]['short_name']}{meta['channel_pos']}"
        sensor_name = short_lookup[meta["sensor_id"]]["long_name"]
        sensor_units = short_lookup[meta["sensor_id"]]["units"]
        coefs = meta["sensor_info"]

        ### Temperature block
        if meta["sensor_id"] == "55":
            log.info(f"Processing Sensor ID: {meta['sensor_id']}, {sensor_name}")
            converted_df[col] = sbe_eq.sbe3(raw_df[meta["column"]], coefs)
            if meta["list_id"] == 0:
                t_array = converted_df[col].astype(float)
                log.info(
                    f"\tPrimary temperature first reading: {t_array[0]} {sensor_units}"
                )

        ### Pressure block
        elif meta["sensor_id"] == "45":
            log.info(f"Processing Sensor ID: {meta['sensor_id']}, {sensor_name}")
            converted_df[col] = sbe_eq.sbe9(raw_df[meta["column"]], t_probe, coefs)
            if meta["list_id"] == 2:
                p_array = converted_df[col].astype(float)
                log.info(f"\tPressure first reading:  {p_array[0]} {sensor_units}")

        ### Conductivity block
        elif meta["sensor_id"] == "3":
            log.info(f"Processing Sensor ID: {meta['sensor_id']}, {sensor_name}")
            converted_df[col] = sbe_eq.sbe4(
                raw_df[meta["column"]], t_array, p_array, coefs
            )
            if meta["list_id"] == 1:
                c_array = converted_df[col].astype(float)
                log.info(f"\tPrimary cond first reading: {c_array[0]} {sensor_units}")

        ### Oxygen block
        elif meta["sensor_id"] == "38":
            log.info(f"Processing Sensor ID: {meta['sensor_id']}, {sensor_name}")
            V_corrected = oxy_eq.sbe43_hysteresis_voltage(
                raw_df[meta["column"]], p_array, coefs
            )
            converted_df[col] = oxy_eq.sbe43(
                V_corrected,
                p_array,
                t_array,
                c_array,
                coefs,
                lat=meta_df["GPSLAT"],
                lon=meta_df["GPSLON"],
            )
            converted_df["CTDOXYVOLTS"] = raw_df[meta["column"]]

        ### Fluorometer Seapoint block
        elif meta["sensor_id"] == "11":
            log.info(f"Processing Sensor ID: {meta['sensor_id']}, {sensor_name}")
            converted_df[col] = aux_eq.seapoint_fluoro(raw_df[meta["column"]], coefs)

        ### Salinity block
        elif meta["sensor_id"] == "1000":
            log.info(f"Processing Sensor ID: {meta['sensor_id']}, {sensor_name}")
            converted_df[col] = gsw.SP_from_C(c_array, t_array, p_array)

        ### Altimeter block
        elif meta["sensor_id"] == "0":
            log.info(f"Processing Sensor ID: {meta['sensor_id']}, {sensor_name}")
            converted_df[col] = aux_eq.sbe_altimeter(raw_df[meta["column"]], coefs)

        ### Rinko block
        elif meta["sensor_id"] == "61":
            if meta["sensor_info"]["SensorName"] in ("RinkoO2V", "RINKO", "RINKOO2", "Rinko02"):
                log.info("Processing Rinko O2")
                # hysteresis correct then pass through voltage (see Uchida, 2010)
                coefs = {"H1": 0.0065, "H2": 5000, "H3": 2000, "offset": 0}
                converted_df[col] = oxy_eq.sbe43_hysteresis_voltage(
                    raw_df[meta["column"]],
                    p_array,
                    coefs,
                )
            elif meta["sensor_info"]["SensorName"] in ("RinkoT"):
                log.info("Processing Rinko T")
                converted_df[col] = raw_df[meta["column"]]

        ### Aux block
        else:
            log.info(f"Passing along Sensor ID: {meta['sensor_id']}, {sensor_name}")
            converted_df[col] = raw_df[meta["column"]]

    # Set the column name for the index
    converted_df.index.name = "index"

    log.info("Joining metadata dataframe with converted data...")
    converted_df = converted_df.join(meta_df)
    log.info("Success!")

    # return the converted data as a dataframe
    return converted_df

def to_temperature(raw, manufacturer, sensor, coefs):
    """
    Wrapper to convert raw temperature output to scientific units using appropriate
    conversion function.

    Parameters
    ----------
    raw : array-like
        Raw temperature output
    manufacturer : str
        Name of sensor manufacturer (e.g. "seabird", "RBR")
    sensor : str
        Name of sensor (e.g. for SBE "3", "4", "9")
    coefs : dict
        Dictionary of calibration coefficient from cal sheet

    Returns
    -------
    temperature : array-like
        Converted temperature (ITS-90)

    """
    # potential example of flexible/user-friendly conversion function

    if manufacturer.lower() in ["sbe", "seabird", "sea-bird"]:
        if sensor == "3":
            # sbe_eq.sbe3()
            pass
        elif sensor == "35":
            # ?
            pass

    elif manufacturer.lower() == "rbr":
        # will need to create equations_rbr module
        if sensor.lower() == "duo":
            pass
        elif sensor.lower() == "concerto":
            pass
def _trim_soak_period(df=None):
    """
    1) Find pump on/off patterns
    2) Select pump_on=True group with largest pressure recording
    3) Find soak period before start of downcast
    4) Trim cast, return everything after top of cast (i.e. minimum pressure)
    """
    df_list = [
        g for i, g in df.groupby(df["pump_on"].ne(df["pump_on"].shift()).cumsum())
    ]
    df_pump_on_list = [df for df in df_list if df["pump_on"].all()]
    df_cast = df_pump_on_list[np.argmax([df["CTDPRS"].max() for df in df_pump_on_list])]
    df_cast = df_cast.reset_index(drop=True)
    # next fn deals w/ edge cases, leave as is for now
    df_cast = _find_last_soak_period(df_cast)
    start_ind = df_cast.loc[: len(df) // 4, "CTDPRS"].argmin()
    df_trimmed = df_cast[start_ind:].reset_index(drop=True).copy()

    return df_trimmed


def _find_last_soak_period(df_cast, time_bin=8, P_surface=2, P_downcast=50):
    """
    Find the soak period before the downcast starts.

    The algorithm is tuned for repeat hydrography work, specifically US GO-SHIP
    parameters. This assumes the soak depth will be somewhere between 10 and 30
    meters, the package will sit at the soak depth for at least 20 to 30 seconds
    before starting ascent to the surface and descent to target depth.

    The algorithm is not guaranteed to catch the exact start of the soak period,
    but within a minimum period of time_bin seconds(?) from end of the soak if
    the soak period assumption is valid. This should be shorter than the total
    soak period time, and able to catch the following rise and descent of the
    package that signals the start of the cast.

    The algorithm has been designed to handle four general cases of casts:
        * A routine cast with pumps turning on in water and normal soak
        * A cast where the pumps turn on in air/on deck
        * A cast where the pumps turn on and off due to rosette coming out of water
        * A cast where there are multiple stops on the downcast to the target depth

    Parameters
    ----------
    df_cast : DataFrame
        DataFrame of the entire cast, from deckbox on to deckbox off
    time_bin : integer, optional
        Number of seconds to bin average for descent rate calculation
    P_surface : integer, optional
        Minimum surface pressure threshold required to look for soak depth
        (2 dbar was chosen as an average rosette is roughly 1.5 to 2 meters tall)
    P_downcast : integer, optional
        Minimum pressure threshold required to assume downcast has started
        (50 dbar has been chosen as double the deep soak depth of 20-30 dbar)

    Returns
    -------
    df_cast_trimmed : DataFrame
        DataFrame starting within time_bin seconds of the last soak period.
    """
    # Validate user input
    if time_bin <= 0:
        raise ValueError("Time bin value should be positive whole seconds.")
    if P_downcast <= 0:
        raise ValueError(
            "Starting downcast pressure threshold must be positive integers."
        )
    if P_downcast < P_surface:
        raise ValueError(
            "Starting downcast pressure threshold must be greater \
                        than surface pressure threshold."
        )

    # If pumps have not turned on until in water, return DataFrame
    if df_cast.iloc[0]["CTDPRS"] > P_surface:
        return df_cast

    # Bin the data by time, and compute the average rate of descent
    df_cast["index"] = df_cast.index  # needed at end to identify start_idx
    df_cast["bin"] = pd.cut(
        df_cast.index,
        np.arange(df_cast.index[0], df_cast.index[-1], time_bin * 24),
        labels=False,
        include_lowest=True,
    )
    df_binned = df_cast.groupby("bin").mean()

    # Compute difference of descent rates and label bins
    df_binned["dP"] = df_binned["CTDPRS"].diff().fillna(0).round(0)
    df_binned["movement"] = pd.cut(
        df_binned["dP"], [-1000, -0.5, 0.5, 1000], labels=["up", "stop", "down"]
    )

    # Find all periods where the rosette is not moving
    df_group = df_binned.groupby(
        df_binned["movement"].ne(df_binned["movement"].shift()).cumsum()
    )
    df_list = [g for i, g in df_group]

    # Find last soak period before starting descent to target depth
    def find_last(df_list, P_downcast):
        for idx, df in enumerate(df_list):
            if df["CTDPRS"].max() < P_downcast:
                # make sure it's soak, not a stop to switch to autocast (i.e. A20 2021)
                if df.max()["movement"] == "stop" and len(df) > 1:
                    last_idx = idx
            else:
                return last_idx
        return last_idx

    # Trim off everything before last soak
    start_idx = int(df_list[find_last(df_list, P_downcast)].head(1)["index"])
    df_cast_trimmed = df_cast.loc[start_idx:].reset_index()

    return df_cast_trimmed


def ctd_align(inMat=None, col=None, time=0.0):
    """ctd_align function

    Function takes full NUMPY ndarray with predefined dtype array
    and adjusts time of sensor responce and water flow relative to
    the time frame of temperature sensor.

    Originally written by Courtney Schatzman, docstring by Joseph Gum.
    Need to generate alignment plots in order to properly use ctd_align.

    Args:
        param1 (ndarray): inMat, numpy ndarray with dtype array
        param2 (float): col, column to apply time advance to.
        param3 (float): time, advance in seconds to apply to raw data.

    Returns:
        Narray: The return value is ndarray with adjusted time of parameter
          specified.

    """
    # Num of frames per second.
    fl = 24

    if (inMat is not None) & (col is not None) & (time > 0.0):
        # Time to advance
        advnc = int(fl * time)
        tmp = np.arange(advnc, dtype=np.float)
        last = inMat[col][len(inMat) - 1]
        tmp.fill(float(last))
        inMat[col] = np.concatenate((inMat[col][advnc:], tmp))

    return inMat


