import xml.etree.ElementTree as ET
import struct
import sys
import os
import json

import sbe_reader as sbe_rd
import sbe_equations_dict as sbe_eq

DEBUG = False

def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def sbe_xml_reader_a(file):
    """Function to read .XMLCON file from Seabird.

    Input:
    file: .XMLCON file produced from Seasave 7.

    Output:
    Dictionary of dictionary of sensors.

    """
    tree = ET.parse(file)
    root = tree.getroot()

    """Pokedex is a dict of {Sensor index numbers from the config:sensor info}
    Assume that sensor index number is the order the sensors have been entered into the file.
    Therefore, it will be Frequency instruments first, then Voltage instruments.
    Call their index number (starting at 0) in order to pull out the info.

    """
    pokedex = {}
    for x in root.iter('Sensor'):
        """Start creating single sensor dictionary."""
        bulbasaur = {}
        bulbasaur['SensorID'] = x.attrib['SensorID']
        #load all values into dict - beware of overwriting #NEED TO FIX
        for children in x:
            for y in children.iter():
                bulbasaur[y.tag] = float_convert(y.text)
            """Add sensor to big dictionary."""
            pokedex[x.attrib['index']] = bulbasaur
    return pokedex

def float_convert(string):
    try:
        return float(string)
    except:
        return string

def cnv_handler_2(hex_file, xmlcon_file, debug=False):
    """Handler to deal with converting eng. data to sci units automatically.
    When not given a format file/json, default to putting out data in order of instruments in xmlcon.
    Format file makes assumptions: duplicate sensors are ranked by channel they use, with lower value more important.
    Ex: For two SBE 3, the temp. sensor on freq. channel 1 is primary temp, and temp. sensor on channel 4 is secondary.

    After reading in format/XMLCON, determine order to put out data.
    Read sensor dictionary, pull out SensorID number, then match with correct method.

    Read in

    VALUES HARDCODED, TRY TO SETUP A DICT TO MAKE NEATER

    """

    global DEBUG
    DEBUG = debug

    debugPrint("Verifying input files... ", end='')

    if not os.path.isfile(hex_file):
        errPrint("Error: .hex file", hex_file, "not found.")
        return False

    if not os.path.isfile(xmlcon_file):
        errPrint("Error: .XMLCOM file", xmlcon_file, "not found.")
        return False

    debugPrint("Success!")


    debugPrint('Reading/Parsing input files...')
    sbe_reader = sbe_rd.SBEReader.from_paths(hex_file, xmlcon_file)

    debugPrint('Retrieving sensor info from xmlcon file...')
    sensor_info = sbe_xml_reader_a(xmlcon_file)
    #debugPrint('Sensor Info:', json.dumps(sensor_info, indent=2))

    #namesplit = xmlcon_file.split('.')
    #namesplit = namesplit[0] + '.converted'
    #outputFile = os.path.splitext(xmlcon_file)[0] + '.converted'
    #debugPrint('Output will be saved to:', outputFile)

    #sensor_dictionary = sbe_xml_reader_a('GS3601101.XMLCON')

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
    processed_data = []

    #Temporary arrays to hold sci_data in order to compute following sci_data (pressure, cond, temp, etc)
    t_array = []
    p_array = []
    c_array = []
    k_array = []

    #lookup table for sensor data
    ###DOUBLE CHECK TYPE IS CORRECT###
    short_lookup = {
        '55':{'short_name': 't', 'long_name':'SBE 3+ Temperature', 'units': 'C', 'type': 'float64'},
        '45':{'short_name': 'p', 'long_name':'SBE 9+ Pressure', 'units': 'dbar', 'type': 'float64'},
        '3':{'short_name': 'c', 'long_name':'SBE 4 Conductivity', 'units': 'S/m', 'type':'float64'},
        '38':{'short_name': 'o', 'long_name':'SBE 43 Oxygen', 'units': 'ml/l', 'type':'float64'},
        '11':{'short_name': 'fluoro', 'long_name':'Seapoint Fluorometer', 'units': 'ug/l', 'type':'float64'},
        '27':{'short_name': 'empty', 'long_name':'empty', 'units':'NA', 'type':'NA'},
        '0':{'short_name': 'alti', 'long_name':'Altitude', 'units':'m', 'type':'float64'},
        '71':{'short_name': 'cstar', 'long_name':'CStar', 'units': 'ug/l', 'type':'float64'},
        '61':{'short_name': 'u_def', 'long_name':'user defined', 'units':'V', 'type':'float64'},
        '1000':{'short_name': 'sal', 'long_name':'Salinity (C1 T1)', 'units':'PSU', 'type':'float64'}
    }

    ######
    # The following are definitions for every key in the dict below:
    #
    # sensor_id = number assigned by SBE for identification in XML
    # list_id = place in XML array by SBE for determining which sensor is which, alternatively channel number (freq+volt)
    # channel_pos = is it the first, second, third, etc sensor of its type in the data file, aux sensors default to 0
    # ranking = data processing ranking - temp first, then pressure, then conductivity, then oxygen, then aux
    # data = eng units to be converted to sci units
    # sensor_info = xml sensor info to convert from eng units to sci units
    ######

    for i, x in enumerate(sensor_info):
        #print(i, sensor_dictionary[str(i)]['SensorID'])
        sensor_id = sensor_info[str(i)]['SensorID']

        #temp block
        if str(sensor_id) == '55':
            temp_counter += 1
            queue_metadata.append({'sensor_id': '55', 'list_id': i, 'channel_pos': temp_counter, 'ranking': 1, 'data': sbe_reader.parsed_scans[:,i], 'sensor_info':sensor_info[str(i)] })

        #cond block
        elif str(sensor_id) == '3':
            cond_counter += 1
            queue_metadata.append({'sensor_id': '3', 'list_id': i, 'channel_pos': cond_counter, 'ranking': 3, 'data': sbe_reader.parsed_scans[:,i], 'sensor_info':sensor_info[str(i)]})

        #pressure block
        elif str(sensor_id) == '45':
            queue_metadata.append({'sensor_id': '45', 'list_id': i, 'channel_pos': '', 'ranking': 2, 'data': sbe_reader.parsed_scans[:,i], 'sensor_info':sensor_info[str(i)]})

        #oxygen block
        elif str(sensor_id) == '38':
            oxygen_counter += 1
            queue_metadata.append({'sensor_id': '38', 'list_id': i, 'channel_pos': oxygen_counter, 'ranking': 5, 'data': sbe_reader.parsed_scans[:,i], 'sensor_info':sensor_info[str(i)]})

        #aux block
        else:
            queue_metadata.append({'sensor_id': sensor_id, 'list_id': i, 'channel_pos': '', 'ranking': 6, 'data': sbe_reader.parsed_scans[:,i], 'sensor_info':sensor_info[str(i)]})

    #a temporary block in order to append basic salinity (t1, c1) to file. If additional salinity is needed (different combinations), it'll need a full reworking
    queue_metadata.append({'sensor_id': '1000', 'list_id': 1000, 'channel_pos':'', 'ranking': 4, 'data': '', 'sensor_info':''})

    queue_metadata = sorted(queue_metadata, key = lambda sensor: sensor['ranking'])

    #queue sorting forces it to be in order, so we don't worry about order here
    #assumes first channel for each sensor is primary for computing following data, rework to accept file to determine which is primary
    while queue_metadata:
        temp_meta = queue_metadata.pop(0)

        ###Temperature block
        if temp_meta['sensor_id'] == '55':
            debugPrint('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            temp_meta['sci_data'] = sbe_eq.temp_its90_dict(temp_meta['sensor_info'], temp_meta['data'])
            if temp_meta['list_id'] == 0:
                t_array = temp_meta['sci_data']
                k_array = [273.15+celcius for celcius in t_array]
                debugPrint('\tPrimary temperature used:', t_array[0], short_lookup[temp_meta['sensor_id']]['units'])
            processed_data.append(temp_meta)
            #debugPrint('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'], short_lookup[temp_meta['sensor_id']]['long_name'])

        ### Pressure block
        elif temp_meta['sensor_id'] == '45':
            debugPrint('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            temp_meta['sci_data'] = sbe_eq.pressure_dict(temp_meta['sensor_info'], temp_meta['data'], t_array)
            if temp_meta['list_id'] == 2:
                p_array = temp_meta['sci_data']
                debugPrint('\tPressure used:', p_array[0], short_lookup[temp_meta['sensor_id']]['units'])
            processed_data.append(temp_meta)
            #debugPrint('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'], short_lookup[temp_meta['sensor_id']]['long_name'])

        ### Conductivity block
        elif temp_meta['sensor_id'] == '3':
            debugPrint('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            temp_meta['sci_data'] = sbe_eq.cond_dict(temp_meta['sensor_info'], temp_meta['data'], t_array, p_array)
            if temp_meta['list_id'] == 1:
                c_array = temp_meta['sci_data']
                debugPrint('\tPrimary cond used:', c_array[0], short_lookup[temp_meta['sensor_id']]['units'])
            processed_data.append(temp_meta)
            #debugPrint('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'], short_lookup[temp_meta['sensor_id']]['long_name'])

        ### Oxygen block
        elif temp_meta['sensor_id'] == '38':
            debugPrint('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            temp_meta['sci_data'] = sbe_eq.oxy_dict(temp_meta['sensor_info'], p_array, k_array, t_array, c_array, temp_meta['data'])
            processed_data.append(temp_meta)
            #debugPrint('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'], short_lookup[temp_meta['sensor_id']]['long_name'])

        ### Fluorometer Seapoint block
        elif temp_meta['sensor_id'] == '11':
            debugPrint('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            temp_meta['sci_data'] = sbe_eq.fluoro_seapoint_dict(temp_meta['sensor_info'], temp_meta['data'])
            processed_data.append(temp_meta)
            #debugPrint('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'], short_lookup[temp_meta['sensor_id']]['long_name'])

        ###Salinity block
        elif temp_meta['sensor_id'] == '1000':
            debugPrint('Processing Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            temp_meta['sci_data'] = sbe_eq.sp_dict(c_array, t_array, p_array)
            processed_data.append(temp_meta)
            #debugPrint('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'], short_lookup[temp_meta['sensor_id']]['long_name'])

        ### Aux block
        else:
            debugPrint('Skipping Sensor ID:', temp_meta['sensor_id'] + ',', short_lookup[temp_meta['sensor_id']]['long_name'])
            temp_meta['sci_data'] = temp_meta['data']
            processed_data.append(temp_meta)
            #debugPrint('Currently skipping (not processing, raw voltages only) sensor list_id: ', temp_meta['sensor_id'], short_lookup[temp_meta['sensor_id']]['long_name'])

    ### Create a single unified object with all data in there
    header_string = []
    header_1 = ''
    header_2 = ''
    header_3 = ''

    """Start writing a .csv file with extension .converted

    First part - compose the header.
    """

    output = ''

    debugPrint('Compiling header information')
    data_list_of_lists = []
    for x in processed_data:
        header_string.append(x['sensor_id'])
        data_list_of_lists.append(x['sci_data'])
        try:
            header_1 = header_1 + '{0}{1},'.format(short_lookup[x['sensor_id']]['short_name'], x['channel_pos'])
        except:
            header_1 = header_1 + 'error,'
            errPrint('Error in lookup table: channel_pos for sensor ID:', x['sensor_id'])

        try:
            header_2 = header_2 + '{0},'.format(short_lookup[x['sensor_id']]['units'])
        except:
            header_2 = header_2 + 'error,'
            errPrint('Error in lookup table: units for sensor ID:', x['sensor_id'])

        try:
            header_3 = header_3 + '{0},'.format(short_lookup[x['sensor_id']]['type'])
        except:
            header_3 = header_3 + 'error,'
            errPrint('Error in lookup table: type for sensor ID:', x['sensor_id'])

    ##### ------------HACKY DATETIME INSERTION------------ #####
    #assumes date/time will always be at end, and adds header accordingly
    #should be rewritten to have cleaner integration with rest of code
    header_1 = header_1 + 'lat,lon,new_pos,nmea_time,scan_time,bottle_fire'
    header_2 = header_2 + 'dec_deg,dec_deg,boolean,ISO8601,ISO8601,boolean'
    header_3 = header_3 + 'float64,float64,bool_,datetime,datetime,bool_'

    ### pos/time/date block
    data_list_of_lists.append(sbe_reader.parsed_scans[:,(sbe_reader.parsed_scans.shape[1]-1)])
    ##### ----------HACKY DATETIME INSERTION END---------- #####

    debugPrint('Transposing data')
    transposed_data = zip(*data_list_of_lists)

    """Write header and body of .csv"""

#    with open(namesplit, 'w') as f:
#    with open(outputFile, 'w') as f:
#        f.write(header_1.rstrip(',') + '\n')
#        f.write(header_2.rstrip(',') + '\n')
#        f.write(header_3.rstrip(',') + '\n')
#        for x in transposed_data:
#            f.write(','.join([str(y) for y in x]) + '\n')

#    with open(namesplit, 'w') as f:
#    with open(outputFile, 'w') as f:
    output += header_1.rstrip(',') + '\n'
    output += header_2.rstrip(',') + '\n'
    output += header_3.rstrip(',') + '\n'
    for x in transposed_data:
        output += ','.join([str(y) for y in x]) + '\n'

#    print('Done, output saved to:', outputFile)
    debugPrint('Processing complete')
    return output
