import xml.etree.ElementTree as ET
import struct
import sys

import sbe_reader as sbe_rd
import sbe_equations_dict as sbe_eq

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

def cnv_handler_2(hex_file, xmlcon_file):
    """Handler to deal with converting eng. data to sci units automatically.
    When not given a format file/json, default to putting out data in order of instruments in xmlcon.
    Format file makes assumptions: duplicate sensors are ranked by channel they use, with lower value more important.
    Ex: For two SBE 3, the temp. sensor on freq. channel 1 is primary temp, and temp. sensor on channel 4 is secondary.

    After reading in format/XMLCON, determine order to put out data.
    Read sensor dictionary, pull out SensorID number, then match with correct method.

    Read in

    VALUES HARDCODED, TRY TO SETUP A DICT TO MAKE NEATER

    """
    sbe_reader = sbe_rd.SBEReader.from_paths(hex_file, xmlcon_file)

    sensor_info = sbe_xml_reader_a(xmlcon_file)
    namesplit = xmlcon_file.split('.')
    namesplit = namesplit[0] + '.converted'
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
    short_lookup = {
        '55':{'short_name': 't', 'units': 'C', 'type': 'float'},
        '45':{'short_name': 'p', 'units': 'dbar', 'type': 'float'},
        '3':{'short_name': 'c', 'units': 'S/m', 'type':'float'},
        '38':{'short_name': 'o', 'units': 'ml/l', 'type':'float'},
        '11':{'short_name': 'fluoro', 'units': 'ug/l', 'type':'float'},
        '27':{'short_name': 'empty', 'units':'NA', 'type':'NA'},
        '0':{'short_name': 'alti', 'units':'m', 'type':'float'},
        '71':{'short_name': 'cstar', 'units': 'ug/l', 'type':'float'},
        '61':{'short_name': 'u_def', 'units':'V', 'type':'float'}
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
            queue_metadata.append({'sensor_id': '38', 'list_id': i, 'channel_pos': oxygen_counter, 'ranking': 4, 'data': sbe_reader.parsed_scans[:,i], 'sensor_info':sensor_info[str(i)]})

        #aux block
        else:
            queue_metadata.append({'sensor_id': sensor_id, 'list_id': i, 'channel_pos': '', 'ranking': 5, 'data': sbe_reader.parsed_scans[:,i], 'sensor_info':sensor_info[str(i)]})

    queue_metadata = sorted(queue_metadata, key = lambda sensor: sensor['ranking'])

    #queue sorting forces it to be in order, so we don't worry about order here
    #assumes first channel for each sensor is primary for computing following data, rework to accept file to determine which is primary
    while queue_metadata:
        temp_meta = queue_metadata.pop(0)

        ###Temperature block
        if temp_meta['sensor_id'] == '55':
            #print(temp_meta['sensor_info'])
            temp_meta['sci_data'] = sbe_eq.temp_its90_dict(temp_meta['sensor_info'], temp_meta['data'])
            if temp_meta['list_id'] == 0:
                t_array = temp_meta['sci_data']
                k_array = [273.15+celcius for celcius in t_array]
                print('Primary temperature used: ', t_array[0])
            processed_data.append(temp_meta)
            print('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'])

        ### Pressure block
        elif temp_meta['sensor_id'] == '45':
            temp_meta['sci_data'] = sbe_eq.pressure_dict(temp_meta['sensor_info'], temp_meta['data'], t_array)
            if temp_meta['list_id'] == 2:
                p_array = temp_meta['sci_data']
                print('Pressure used: ', p_array[0])
            processed_data.append(temp_meta)
            print('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'])

        ### Conductivity block
        elif temp_meta['sensor_id'] == '3':
            temp_meta['sci_data'] = sbe_eq.cond_dict(temp_meta['sensor_info'], temp_meta['data'], t_array, p_array)
            if temp_meta['list_id'] == 1:
                c_array = temp_meta['sci_data']
                print('Primary cond used: ', c_array[0])
            processed_data.append(temp_meta)
            print('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'])

        ### Oxygen block
        elif temp_meta['sensor_id'] == '38':
            temp_meta['sci_data'] = sbe_eq.oxy_dict(temp_meta['sensor_info'], p_array, k_array, t_array, c_array, temp_meta['data'])
            processed_data.append(temp_meta)
            print('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'])

        ### Fluorometer Seapoint block
        elif temp_meta['sensor_id'] == '11':
            temp_meta['sci_data'] = sbe_eq.fluoro_seapoint_dict(temp_meta['sensor_info'], temp_meta['data'])
            processed_data.append(temp_meta)
            print('Processed ', temp_meta['ranking'], temp_meta['list_id'], temp_meta['sensor_id'])

        ### Aux block
        else:
            temp_meta['sci_data'] = temp_meta['data']
            processed_data.append(temp_meta)
            print('Currently skipping (not processing, raw voltages only) sensor list_id: ', temp_meta['list_id'])

    ### Create a single unified object with all data in there
    header_string = []
    header_1 = ''
    header_2 = ''
    header_3 = ''

    data_list_of_lists = []
    for x in processed_data:
        header_string.append(x['sensor_id'])
        data_list_of_lists.append(x['sci_data'])
        try:
            header_1 = header_1 + '{0}{1},'.format(short_lookup[x['sensor_id']]['short_name'], x['channel_pos'])
            header_2 = header_2 + '{0},'.format(short_lookup[x['sensor_id']]['units'])
            header_3 = header_3 + '{0},'.format(short_lookup[x['sensor_id']]['type'])
        except:
            print(None)

    transposed_data = zip(*data_list_of_lists)

    #print(header_1.rstrip(','))
    #print(header_2)
    #print(header_3)
    #print(header_string)

    """
    #IT WORKS
    for x in transposed_data:
        print(','.join([str(y) for y in x]))
    """

    with open(namesplit, 'w') as f:
        f.write(header_1.rstrip(',') + '\n')
        f.write(header_2.rstrip(',') + '\n')
        f.write(header_3.rstrip(',') + '\n')
        for x in transposed_data:
            f.write(','.join([str(y) for y in x]) + '\n')

    print('Done, look for file with ".converted" extension')
    return None
