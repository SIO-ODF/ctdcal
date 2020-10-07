from pathlib import Path

import config as cfg
import gsw
import pandas as pd

import ctdcal.process_bottle as btl
import ctdcal.process_ctd as process_ctd
import ctdcal.sbe_equations_dict as sbe_eq
import ctdcal.sbe_reader as sbe_rd

DEBUG = False

#lookup table for sensor data
###DOUBLE CHECK TYPE IS CORRECT###
short_lookup = {
    '55':{'short_name': 'CTDTMP', 'long_name':'SBE 3+ Temperature', 'units': 'ITS-90', 'type': 'float64'},
    '45':{'short_name': 'CTDPRS', 'long_name':'SBE 9+ Pressure', 'units': 'DBAR', 'type': 'float64'},
    '3':{'short_name': 'CTDCOND', 'long_name':'SBE 4 Conductivity', 'units': 'MSPCM', 'type':'float64'},
    '38':{'short_name': 'CTDOXY', 'long_name':'SBE 43 Oxygen', 'units': 'MLPL', 'type':'float64'},
    #'38':{'short_name': 'CTDOXYVOLTS', 'long_name':'SBE 43 Oxygen Volts', 'units': '0-5VDC', 'type':'float64'},
    '11':{'short_name': 'FLUOR', 'long_name':'Seapoint Fluorometer', 'units': '0-5VDC', 'type':'float64'},
    '27':{'short_name': 'FREE', 'long_name':'empty', 'units':'NA', 'type':'NA'},
    '0':{'short_name': 'ALT', 'long_name':'Altitude', 'units':'M', 'type':'float64'},
    '71':{'short_name': 'CTDXMISS', 'long_name':'CStar', 'units': '0-5VDC', 'type':'float64'},
    '61':{'short_name': 'U_DEF', 'long_name':'user defined', 'units':'0-5VDC', 'type':'float64'},
    '1000':{'short_name': 'CTDSAL', 'long_name':'Salinity (C1 T1)', 'units':'PSU', 'type':'float64'},
    '20':{'short_name': 'CTDFLUOR', 'long_name':'WetlabECO_AFL_FL_Sensor', 'units':'0-5VDC', 'type':'float64'}, #check short_name later
    '42':{'short_name':'PAR', 'long_name':'PAR/Irradiance, Biospherical/Licor', 'units':'0-5VDC', 'type':'float64'},
    '51':{'short_name':'REF_PAR', 'long_name':'Surface PAR/Irradiance, Biospherical/Licor', 'units':'0-5VDC', 'type':'float64'},
    '70':{'short_name': 'CTDBACKSCATTER', 'long_name': 'WetlabECO_BB_Sensor', 'units':'0-5VDC', 'type':'float64'}
}


def hex_to_ctd(ssscc_list, debug=False):
    # TODO: add (some) error handling from odf_convert_sbe.py
    """
    Convert raw CTD data and export to .pkl files.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to convert
    debug : bool, optional
        Display verbose messages

    Returns
    -------

    """
    # TODO: use logger module instead
    print('Converting .hex files')
    for ssscc in ssscc_list:
        if not Path(cfg.directory["converted"] + ssscc + ".pkl").exists():
            hexFile = cfg.directory["raw"] + ssscc + ".hex"
            xmlconFile = cfg.directory["raw"] + ssscc + ".XMLCON"
            sbeReader = sbe_rd.SBEReader.from_paths(hexFile, xmlconFile)
            converted_df = convertFromSBEReader(sbeReader, debug=debug)
            converted_df.to_pickle(cfg.directory["converted"] + ssscc + ".pkl")

    return True


def make_time_files(ssscc_list):
    print("Generating time.pkl files")
    for ssscc in ssscc_list:
        if not Path(cfg.directory["time"] + ssscc + "_time.pkl").exists():
            converted_df = pd.read_pickle(cfg.directory["converted"] + ssscc + ".pkl")

            # Trim to times when rosette is in water
            trimmed_df = process_ctd.ondeck_pressure_2(
                converted_df,
                ssscc,
                log_file=cfg.directory["logs"] + "ondeck_pressure.csv",
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
                trimmed_df, window="triangle", parameters=cfg.filter_cols,
            )

            # Trim to downcast
            cast_data = process_ctd.cast_details(
                filter_data, ssscc, log_file=cfg.directory["logs"] + "cast_details.csv",
            )

            cast_data.to_pickle(cfg.directory["time"] + ssscc + "_time.pkl")


def make_btl_mean(ssscc_list, debug=False):
    # TODO: add (some) error handling from odf_process_bottle.py
    """
    Create "bottle mean" files from continuous CTD data averaged at the bottle stops.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to convert
    debug : bool, optional
        Display verbose messages

    Returns
    -------
    boolean
        bottle averaging of mean has finished successfully
    """
    print('Generating btl_mean.pkl files')
    for ssscc in ssscc_list:
        if not Path(cfg.directory["bottle"] + ssscc + "_btl_mean.pkl").exists():
            imported_df = pd.read_pickle(cfg.directory["converted"] + ssscc + ".pkl")
            bottle_df = btl.retrieveBottleData(imported_df, debug=debug)
            mean_df = btl.bottle_mean(bottle_df)
            mean_df.to_pickle(cfg.directory["bottle"] + ssscc + "_btl_mean.pkl")

    return True


def convertFromSBEReader(sbeReader, debug=False):
    """Handler to convert engineering data to sci units automatically.
    Takes SBEReader object that is already connected to the .hex and .XMLCON files.
    Optionally takes a boolean debug flag to specify whether or not to display
    verbose messages to stderr
    """

    global DEBUG
    DEBUG = debug

    # Retrieve parsed scans
    rawData = sbeReader.parsed_scans

    # Convert raw data to dataframe
    raw_df = pd.DataFrame(rawData)
    raw_df.index.name = 'index'
    raw_df = raw_df.apply(pd.to_numeric, errors="ignore")

    # Retrieve Config data
    rawConfig = sbeReader.parsed_config()

    # The meta data field needs to be processed seperately and then joined with the converted_df
    print("Building meta data dataframe... ")
    metaArray = [line.split(',') for line in sbeReader._parse_scans_meta().tolist()]
    metaArrayheaders = sbeReader._breakdown_header()
    meta_df = pd.DataFrame(metaArray)
    #print(meta_df)
    #print(metaArrayheaders[0])
    meta_df.columns = metaArrayheaders[0]
    meta_df.index.name = 'index'

    for i, x in enumerate(metaArrayheaders[0]):
        if not metaArrayheaders[1][i] == 'bool_':
            meta_df[metaArrayheaders[0][i]] = meta_df[metaArrayheaders[0][i]].astype(metaArrayheaders[1][i])
        else:
            meta_df[metaArrayheaders[0][i]] = meta_df[metaArrayheaders[0][i]].str.match('True', na=False)

    print('Success!')

    pressure_temp = meta_df['pressure_temp_int'].tolist()
    #needs to search sensor dictionary, and compute in order:
    #temp, pressure, cond, salinity, oxygen, all aux.
    #run one loop that builds a queue to determine order of processing, must track which column to pull
    #process queue, store results in seperate arrays for reuse later
    #once queue is empty, attach results together according to format order or xmlcon order - structure to keep track
    queue_metadata = []
    results = {}
    temp_counter = 0
    cond_counter = 0
    oxygen_counter = 0
    u_def_counter = 0
    empty_counter = 0
    processed_data = []

    #Temporary arrays to hold sci_data in order to compute following sci_data (pressure, cond, temp, etc)
    t_array = []
    p_array = []
    c_array = []
    k_array = []

    ######
    # The following are definitions for every key in the dict below:
    #
    # sensor_id = number assigned by SBE for identification in XML
    # list_id = place in XML array by SBE for determining which sensor is which, alternatively channel number (freq+volt)
    # channel_pos = is it the first, second, third, etc sensor of its type in the data file, aux sensors default to 0
    # ranking = data processing ranking - temp first, then pressure, then conductivity, then oxygen, then aux
    # column = column in the raw_df containing the engineering units to be converted to sci units
    # sensor_info = xml sensor info to convert from eng units to sci units
    ######

    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']

        #temp block
        if sensor_id == '55':
            temp_counter += 1
            queue_metadata.append({'sensor_id': '55', 'list_id': i, 'channel_pos': temp_counter, 'ranking': 1, 'column': i, 'sensor_info':rawConfig['Sensors'][i]})

        #cond block
        elif str(sensor_id) == '3':
            cond_counter += 1
            queue_metadata.append({'sensor_id': '3', 'list_id': i, 'channel_pos': cond_counter, 'ranking': 3, 'column': i, 'sensor_info':rawConfig['Sensors'][i]})

        #pressure block
        elif str(sensor_id) == '45':
            queue_metadata.append({'sensor_id': '45', 'list_id': i, 'channel_pos': '', 'ranking': 2, 'column': i, 'sensor_info':rawConfig['Sensors'][i]})

        #oxygen block
        elif str(sensor_id) == '38':
            oxygen_counter += 1
            queue_metadata.append({'sensor_id': '38', 'list_id': i, 'channel_pos': oxygen_counter, 'ranking': 5, 'column': i, 'sensor_info':rawConfig['Sensors'][i]})

        #empty block
        elif str(sensor_id) == '27':
            empty_counter += 1
            queue_metadata.append({'sensor_id': '27', 'list_id': i, 'channel_pos': empty_counter, 'ranking': 6, 'column': i, 'sensor_info':rawConfig['Sensors'][i]})

        #u_def block
        elif str(sensor_id) == '61':
            u_def_counter += 1
            queue_metadata.append({'sensor_id': '61', 'list_id': i, 'channel_pos': u_def_counter, 'ranking': 6, 'column': i, 'sensor_info':rawConfig['Sensors'][i]})

        #aux block
        else:
            queue_metadata.append({'sensor_id': sensor_id, 'list_id': i, 'channel_pos': '', 'ranking': 7, 'column': i, 'sensor_info':rawConfig['Sensors'][i]})

    #a temporary block in order to append basic salinity (t1, c1) to file. If additional salinity is needed (different combinations), it'll need a full reworking
    queue_metadata.append({'sensor_id': '1000', 'list_id': 1000, 'channel_pos':'', 'ranking': 4, 'column': '', 'sensor_info':''})

    #queue sorting forces it to be in order, so we don't worry about order here
    #assumes first channel for each sensor is primary for computing following data, rework to accept file to determine which is primary
    queue_metadata = sorted(queue_metadata, key = lambda sensor: sensor['ranking'])

    #empty converted dataframs
    converted_df = pd.DataFrame()

    for temp_meta in queue_metadata:

        column_name = '{0}{1}'.format(short_lookup[temp_meta['sensor_id']]['short_name'], temp_meta['channel_pos'])

        ###Temperature block
        if temp_meta['sensor_id'] == '55':
            print('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = sbe_eq.temp_its90_dict(temp_meta['sensor_info'], raw_df[temp_meta['column']])
            if temp_meta['list_id'] == 0:
                t_array = converted_df[column_name].astype(float)
                k_array = [273.15+celcius for celcius in t_array]
                print('\tPrimary temperature first reading:', t_array[0], short_lookup[temp_meta['sensor_id']]['units'])
            #processed_data.append(temp_meta)

        ### Pressure block
        elif temp_meta['sensor_id'] == '45':
            print('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = sbe_eq.pressure_dict(temp_meta['sensor_info'], raw_df[temp_meta['column']], pressure_temp)
            if temp_meta['list_id'] == 2:
                p_array = converted_df[column_name].astype(float)
                print('\tPressure first reading:', p_array[0], short_lookup[temp_meta['sensor_id']]['units'])
            #processed_data.append(temp_meta)

        ### Conductivity block
        elif temp_meta['sensor_id'] == '3':
            print('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = sbe_eq.cond_dict(temp_meta['sensor_info'], raw_df[temp_meta['column']], t_array, p_array)
            if temp_meta['list_id'] == 1:
                c_array = converted_df[column_name].astype(float)
                print('\tPrimary cond first reading:', c_array[0], short_lookup[temp_meta['sensor_id']]['units'])
            #processed_data.append(temp_meta)

        ### Oxygen block
        elif temp_meta['sensor_id'] == '38':
            print('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = sbe_eq.oxy_dict(temp_meta['sensor_info'], p_array, k_array, t_array, c_array, raw_df[temp_meta['column']])
            converted_df['CTDOXYVOLTS'] = raw_df[temp_meta['column']]
            #processed_data.append(temp_meta)

        ### Fluorometer Seapoint block
        elif temp_meta['sensor_id'] == '11':
            print('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = sbe_eq.fluoro_seapoint_dict(temp_meta['sensor_info'], raw_df[temp_meta['column']])
            #processed_data.append(temp_meta)

        ###Salinity block
        elif temp_meta['sensor_id'] == '1000':
            print('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = gsw.SP_from_C(c_array, t_array, p_array)
            #processed_data.append(temp_meta)

        ###Altimeter block
        elif temp_meta['sensor_id'] == '0':
            print('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = sbe_eq.altimeter_voltage(temp_meta['sensor_info'], raw_df[temp_meta['column']])

        ### Aux block
        else:
            print('Passing along Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            converted_df[column_name] = raw_df[temp_meta['column']]
            #processed_data.append(temp_meta)

    # Set the column name for the index
    converted_df.index.name = 'index'

    print("Joining meta data dataframe with converted data... ")
    converted_df = converted_df.join(meta_df)
    print('Success!')

    # return the converted data as a dataframe
    return converted_df
