from pathlib import Path

import gsw
import logging
import numpy as np
import pandas as pd

from . import get_ctdcal_config
from . import equations_sbe as sbe_eq
from . import process_bottle as btl
from . import process_ctd as process_ctd
from . import sbe_reader as sbe_rd

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

# TODO: move this to a separate file?
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
}


def hex_to_ctd(ssscc_list):
    # TODO: add (some) error handling from odf_convert_sbe.py
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
        if not Path(cfg.dirs["converted"] + ssscc + ".pkl").exists():
            hexFile = cfg.dirs["raw"] + ssscc + ".hex"
            xmlconFile = cfg.dirs["raw"] + ssscc + ".XMLCON"
            sbeReader = sbe_rd.SBEReader.from_paths(hexFile, xmlconFile)
            converted_df = convertFromSBEReader(sbeReader)
            converted_df.to_pickle(cfg.dirs["converted"] + ssscc + ".pkl")

    return True


def make_time_files(ssscc_list):
    log.info("Generating time.pkl files")
    for ssscc in ssscc_list:
        if not Path(cfg.dirs["time"] + ssscc + "_time.pkl").exists():
            converted_df = pd.read_pickle(cfg.dirs["converted"] + ssscc + ".pkl")

            # Remove any pressure spikes
            bad_rows = converted_df["CTDPRS"].abs() > 6500
            if bad_rows.any():
                log.debug(f"{ssscc}: {bad_rows.sum()} bad pressure points removed.")
            converted_df.loc[bad_rows, :] = np.nan
            converted_df.interpolate(limit=24, limit_area="inside", inplace=True)

            # Trim to times when rosette is in water
            trimmed_df = process_ctd.remove_on_deck(
                converted_df,
                ssscc,
                log_file=cfg.dirs["logs"] + "ondeck_pressure.csv",
            )

            # # TODO: switch to loop instead, e.g.:
            # align_cols = [cfg.column[x] for x in ["c1", "c2"]]  # "dopl" -> "CTDOXY1"

            # if not c1_col in raw_data.dtype.names:
            #     print('c1_col data not found, skipping')
            # else:
            #     raw_data = process_ctd.ctd_align(raw_data, c1_col, float(tc1_align))
            # if not c2_col in raw_data.dtype.names:
            #     print('c2_col data not found, skipping')
            # else:
            #     raw_data = process_ctd.ctd_align(raw_data, c2_col, float(tc2_align))
            # if not dopl_col in raw_data.dtype.names:
            #     print('do_col data not found, skipping')
            # else:
            #     raw_data = process_ctd.ctd_align(raw_data, dopl_col, float(do_align))

            # TODO: add despike/wild edit filter (optional?)

            # Filter data
            filter_data = process_ctd.raw_ctd_filter(
                trimmed_df,
                window="triangle",
                parameters=cfg.filter_cols,
            )

            # Trim to downcast
            cast_data = process_ctd.cast_details(
                filter_data,
                ssscc,
                log_file=cfg.dirs["logs"] + "cast_details.csv",
            )

            cast_data.to_pickle(cfg.dirs["time"] + ssscc + "_time.pkl")


def make_btl_mean(ssscc_list):
    # TODO: add (some) error handling from odf_process_bottle.py
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
            bottle_df = btl.retrieveBottleData(imported_df)
            mean_df = btl.bottle_mean(bottle_df)

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


def convertFromSBEReader(sbeReader):
    """Handler to convert engineering data to sci units automatically.
    Takes SBEReader object that is already connected to the .hex and .XMLCON files.
    """

    # Retrieve parsed scans and convert to dataframe
    rawData = sbeReader.parsed_scans
    raw_df = pd.DataFrame(rawData).apply(pd.to_numeric, errors="ignore")
    raw_df.index.name = "index"

    # Metadata needs to be processed seperately and then joined with the converted data
    print("Building metadata dataframe...")
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

    print("Success!")

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
    # TODO: rework to use config.py file to determine which is primary
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
            print("Processing Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            converted_df[col] = sbe_eq.sbe3(raw_df[meta["column"]], coefs)
            if meta["list_id"] == 0:
                t_array = converted_df[col].astype(float)
                print("\tPrimary temperature first reading:", t_array[0], sensor_units)

        ### Pressure block
        elif meta["sensor_id"] == "45":
            print("Processing Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            converted_df[col] = sbe_eq.sbe9(raw_df[meta["column"]], t_probe, coefs)
            if meta["list_id"] == 2:
                p_array = converted_df[col].astype(float)
                print("\tPressure first reading:", p_array[0], sensor_units)

        ### Conductivity block
        elif meta["sensor_id"] == "3":
            print("Processing Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            converted_df[col] = sbe_eq.sbe4(
                raw_df[meta["column"]], t_array, p_array, coefs
            )
            if meta["list_id"] == 1:
                c_array = converted_df[col].astype(float)
                print("\tPrimary cond first reading:", c_array[0], sensor_units)

        ### Oxygen block
        elif meta["sensor_id"] == "38":
            print("Processing Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            # TODO: put some kind of user-enabled flag in config.py, e.g.
            # if cfg["correct_oxy_hysteresis"]:
            V_corrected = sbe_eq.sbe43_hysteresis_voltage(
                raw_df[meta["column"]], p_array, coefs
            )
            converted_df[col] = sbe_eq.sbe43(
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
            print("Processing Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            converted_df[col] = sbe_eq.seapoint_fluoro(raw_df[meta["column"]], coefs)

        ### Salinity block
        elif meta["sensor_id"] == "1000":
            print("Processing Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            converted_df[col] = gsw.SP_from_C(c_array, t_array, p_array)

        ### Altimeter block
        elif meta["sensor_id"] == "0":
            print("Processing Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            converted_df[col] = sbe_eq.sbe_altimeter(raw_df[meta["column"]], coefs)

        elif meta["sensor_id"] == "61":
            if meta["sensor_info"]["SensorName"] == "RinkoO2V":
                print("Processing Rinko O2")
                # hysteresis correct then pass through voltage (see Uchida, 2010)
                coefs = {"H1": 0.0065, "H2": 5000, "H3": 2000, "offset": 0}
                converted_df[col] = sbe_eq.sbe43_hysteresis_voltage(
                    raw_df[meta["column"]],
                    p_array,
                    coefs,
                )

        ### Aux block
        else:
            print("Passing along Sensor ID:", meta["sensor_id"] + ",", sensor_name)
            converted_df[col] = raw_df[meta["column"]]

    # Set the column name for the index
    converted_df.index.name = "index"

    print("Joining metadata dataframe with converted data...")
    converted_df = converted_df.join(meta_df)
    print("Success!")

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


def CR_to_cond(cr, bath_t, ref_t, btl_p):
    """
    Convert AutoSal double conductivity ratio (CR) to conductivity using
    GSW conversion routines.

    Parameters
    ----------
    cr : array-like
        Double conductivity ratio from AutoSal, unitless
    bath_t : array-like
        AutoSal water bath temperature, degrees C
    ref_t : array-like
        CTD temperature at bottle stop, degrees C
    btl_p : array-like
        CTD pressure at bottle stop, dbar

    Returns
    -------
    cond : array-like
        Converted reference conductivity, mS/cm

    """

    salinity = gsw.SP_salinometer((cr / 2.0), bath_t)
    cond = gsw.C_from_SP(salinity, ref_t, btl_p)

    # ignore RunTimeWarning from (np.nan <= 1)
    with np.errstate(invalid="ignore"):
        cond[cond <= 1] = np.nan

    return cond
