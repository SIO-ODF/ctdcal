#!/bin/bash/env
'''
This module is used to parse the xml config and hex file provided by user. It
takes the xml config and extracts the useful settings used during a CTD cast,
and uses that to check and parse the hex file.
'''

import xml.etree.cElementTree as ET
import struct
import re
import sys
import datetime
from pytz import timezone
import numpy as np
import pandas as pd


class SBEReader():
    '''
    Class: SBEReader()
    Description: SBEReader is designed to take the hex and xmlconf files provided
    by the user and parse through both to extract neeeded data for other classes
    that call upon this class.
    Parameter(s): none
    Return(s): none
    '''
    def __init__(self, raw_hex, xml_config):
        '''
        Function: _init_
        Description: Creates a new instance of the SBEReader object using the
        provided raw hex file and xml config file to create the instance.
        Parameter(s): self, raw_hex, xml_config
        Return(s): none
        '''
        self.raw_hex = raw_hex
        self.xml_config = xml_config
        self._parse_config()
        self._load_hex()
        self._check_scan_lengths()

    def _parse_config(self):
        '''
        Function: _parse_config_
        Description: This function takes the xml config and parses the relevant
        information from it.
        Parameter(s): self
        Return(s): none
        '''
        added_checklist = [
            "SurfaceParVoltageAdded",
            "NmeaPositionDataAdded",
            "NmeaDepthDataAdded",
            "NmeaTimeAdded",
            "ScanTimeAdded"
        ]
        suppressed_checklist = [
            "FrequencyChannelsSuppressed",
            "VoltageWordsSuppressed"
        ]
        sensors = {}

        config = ET.fromstring(self.xml_config)
        self.config = {}

        for key in suppressed_checklist:
            try:
                self.config[key] = \
                    int(config.find("./Instrument/{}".format(key)).text)
            except AttributeError:
                raise AttributeError("Could not find {} in \
                    XMLCONF".format(key))
            except ValueError:
                raise ValueError("{} Value is not accurate".format(key))

        for key in added_checklist:
            try:
                self.config[key] = \
                    bool(int(config.find("./Instrument/{}".format(key)).text))
            except AttributeError:
                raise AttributeError("Could not find {} in\
                        XMLCONF".format(key))
            except ValueError:
                raise ValueError("{} Value is not a number".format(key))

        for x in config.iter("Sensor"):
            sensor_dict = {}
            sensor_dict["SensorID"] = x.attrib["SensorID"]

            for children in x:
                for y in children.iter():
                    try:
                        sensor_dict[y.tag] = float(y.text)
                    except:
                        sensor_dict[y.tag] = str(y.text).replace(\
                            "\n", "").replace(" ", "")

                sensors[int(x.attrib["index"])] = sensor_dict

        self.config["Sensors"] = sensors





    def _parse_scans(self):
        '''
        Function: _parse_scans_
        Description: This function takes the relevant information extracted
        from the xml config file and the hex file to get the data collected
        in the hex file produced by the SBE instrument. Returns the data as
        the measurements variable.
        Parameter(s): self
        Return(s): measurements
        '''
        num_frequencies = 5 - self.config["FrequencyChannelsSuppressed"]
        num_voltages = 8 - self.config["VoltageWordsSuppressed"]
        flag_spar = 0

        if self.config["SurfaceParVoltageAdded"]:
            flag_spar = 1

        the_bytes = b"".join(self.raw_bytes)

        unpack_str = "6s" * num_frequencies + "3s" * num_voltages + "2s" *\
            flag_spar + "4s" * flag_spar + "{}s".format(self.scan_length -
                                                        num_voltages * 3 -
                                                        num_frequencies * 6 -
                                                        flag_spar * 6)

        measurements = pd.DataFrame([[int(x, 16) for x in line] for line in\
            struct.iter_unpack(unpack_str, the_bytes)])
        measurements.loc[:, :num_frequencies] = \
            measurements.loc[:, :num_frequencies] / 256
        measurements.loc[:, num_frequencies:num_frequencies + \
                        num_voltages] = (5 * (1 -
                                              (measurements.loc[:,
                                                                num_frequencies
                                                                :num_frequencies
                                                                + num_voltages]
                                               / 4095)))
        return measurements

    def _parse_scans_meta(self):
        '''
        Function: _parse_scans_meta
        Description: This function takes the relevant information extracted
        from the xml config file and the hex file to get the meta data
        collected in the hex file produced by the SBE instrument. Returns
        either measurements_2 or measurements_3 depending on the set scan time
        flag variable. Both return variables contain some form of the meta data
        from the SBE instrumenti after being unpacked by _breakdown().
        Parameter(s): self
        Return(s): measurements_2 or measurements_3
        '''

        num_frequencies = 5 - self.config["FrequencyChannelsSuppressed"]
        num_voltages = 8 - self.config["VoltageWordsSuppressed"]

        string_order = ["scan"]
        flag_nmea_pos = 0
        flag_nmea_depth = 0
        flag_nmea_time = 0
        flag_scan_time = 0
        flag_spar = 0
        flag_pressure_temp = 1
        flag_ctd_status = 1
        flag_modulo = 1

        if self.config["SurfaceParVoltageAdded"]:
            flag_spar = 1
        if self.config["NmeaPositionDataAdded"]:
            flag_nmea_pos = 1
            string_order.append("nmea_pos")
        if self.config["NmeaDepthDataAdded"]:
            flag_nmea_depth = 1
            string_order.append("nmea_depth")
        if self.config["NmeaTimeAdded"]:
            flag_nmea_time = 1
            string_order.append("nmea_time")
        string_order.append("pressure_temp")
        string_order.append("flag_ctd_status")
        string_order.append("modulo")
        if self.config["ScanTimeAdded"]:
            flag_scan_time = 1
            string_order.append("scan_time")

        the_bytes = b"".join(self.raw_bytes)

        unpack_str = (str(num_frequencies * 6 + num_voltages * 3 + flag_spar *
                          6) + "s" +
                      "14s" * flag_nmea_pos +
                      "6s" * flag_nmea_depth +
                      "8s" * flag_nmea_time +
                      "3s" * flag_pressure_temp +
                      "1s" * flag_ctd_status +
                      "2s" * flag_modulo +
                      "8s" * flag_scan_time)

        measurements = ([line for line in struct.iter_unpack(unpack_str,
                                                             the_bytes)])
        measurements_2 = [self._breakdown(line, string_order) for line in\
            measurements]


        if flag_scan_time == 0:
            measurements_3 = self._sbe_time_seq(measurements_2)
            return np.array(measurements_3)

        return np.array(measurements_2)



# TODO: Combine parse_scans and parse_scans_meta functions

#    def _parse_scans(self):
#
#        frequency_offset = 5
#        voltage_offset = 8
#        string_order = ["scan"]
#        # flags to determine how many bytes to break by
#        flag_nmea_pos = 0
#        flag_nmea_depth = 0
#        flag_nmea_time = 0
#        flag_scan_time = 0
#        flag_spar = 0
#        # put in to make code readable later on
#        flag_pressure_temp = 1
#        flag_ctd_status = 1
#        flag_modulo = 1
#
#        # parse the frequencies
#        num_frequencies = frequency_offset -\
#            self.config["FrequencyChannelsSuppressed"]
#        num_voltages = voltage_offset - self.config["VoltageWordsSuppressed"]
#
#        if self.config["SurfaceParVoltageAdded"]:
#            flag_spar = 1
#        if self.config["NmeaPositionDataAdded"]:
#            flag_nmea_pos = 1
#            string_order.append("nmea_pos")
#        # Depth is here for completeness but not implemented, after email chain
#        # showed SBE no longer knows how they did it.
#        if self.config["NmeaDepthDataAdded"]:
#            flag_nmea_depth = 1
#            string_order.append("nmea_depth")
#        if self.config["NmeaTimeAdded"]:
#            flag_nmea_time = 1
#            string_order.append("nmea_time")
#        string_order.append("pressure_temp")
#        string_order.append("flag_ctd_status")
#        string_order.append("modulo")
#        if self.config["ScanTimeAdded"]:
#            flag_scan_time = 1
#            string_order.append("scan_time")
#
#        the_bytes = b"".join(self.raw_bytes)
#        # specify how each line of raw data should be broken down
#        unpack_str = "6s" * num_frequencies +\
#            "3s" * num_voltages +\
#            "2s" * flag_spar +\
#            "4s" * flag_spar +\
#            "{}s".format(self.scan_length - num_voltages * 3 -
#                num_frequencies * 6 - flag_spar * 6)
#        measurements = pd.DataFrame([[\
#            int(x, 16) for x in line] for line in struct.iter_unpack(\
#                unpack_str, the_bytes)])
#
#        # this format is fed directly into _breakdown, so it needs to be
#        # tracked there.
#        meta_unpack_str = (str(\
#            num_frequencies * 6 + num_voltages * 3 + flag_spar * 6) + "s" +\
#            "14s" * flag_nmea_pos +
#            "6s" * flag_nmea_depth +
#            "8s" * flag_nmea_time +
#            "3s" * flag_pressure_temp +
#            "1s" * flag_ctd_status +
#            "2s" * flag_modulo +
#            "8s" * flag_scan_time)
#
#        # this uses pandas magic to just "apply" the needed operations to the
#        # correct part of the array all in one go
#        measurements.loc[:,:num_frequencies] /= 256
#        measurements.loc[:,num_frequencies:num_frequencies + num_voltages] = \
#            (5 * (1 - (measurements.loc[:,num_frequencies:num_frequencies +
#            num_voltages] / 4095)))
#
#        # breaks down line according to columns specified by meta_unpack_str
#        # and converts from hex to decimal MSB first needs to be adjusted for
#        # LSB fields (NMEA time, scan time), break into two lines? move int(x,
#        # 16) outside
#        meta_measurements = [line for line in struct.iter_unpack(
#            meta_unpack_str,the_bytes)]
#        meta_measurements_2 = [self._breakdown(line, string_order) for line in measurements]
#
#        if flag_scan_time == 0:
#            meta_measurements_3 = self._sbe_time_seq(meta_measurements_2)
#            meta_measurements_3 = pd.DataFrame(meta_measurements_3)
#            return meta_measurements_3
#        else:
#            return pd.DataFrame(meta_measurements_2)
#
#        return measurements



    def _breakdown(self, line, ordering):
        '''
        Function: _breakdown
        Description: This function breaks down the meta data collected in
        _parse_scans_meta and 'unpacks' it. The logic on unpacking is set by
        the SBE instrument and has a structure to it.
        Parameter(s): self, line, ordering
        Return(s): output
        '''
        output = ""
        for x, y in zip(ordering, line):
            if x == "nmea_pos":
                tokens = []
                for t in struct.iter_unpack("2s2s2s2s2s2s2s", y):
                    for tt in t:
                        tokens.append(tt)
                output = output + str(self._location_fix(tokens[0], tokens[1],
                                                         tokens[2], tokens[3],
                                                         tokens[4], tokens[5],
                                                         tokens[6])) + ","
            elif x == "nmea_time":
                output = output + str(self._sbe_time(self._reverse_bytes(y),
                                                     "nmea")) + ","
            elif x == "scan_time":
                output = output + str(self._sbe_time(self._reverse_bytes(y),
                                                     "scan"))
            elif x == "flag_ctd_status":
                # make final comma dynamic later on so less errors down the
                # line
                output = output + str(self._pump_status(y)) + "," +\
                    str(self._bottle_fire(y)) + ","
            elif x == "pressure_temp":
                output = output + str(int(y, 16)) + ","
        return output

    def _breakdown_header(self):
        '''
        Function: _breakdown_header
        Description: This function breaks down the header of the xml config and
        extracts relevant information and configure it for use.
        Parameter(s): self
        Return(s): output
        '''

        output = [[], []]

        if self.config["NmeaPositionDataAdded"]:
            output[0] = ["GPSLAT", "GPSLON", "new_fix"]
            output[1] = ["float64", "float64", "bool_"]
        if self.config["NmeaTimeAdded"]:
            output[0].append("nmea_datetime")
            output[1].append("float64")

        output[0].append("pressure_temp_int")
        output[1].append("int_")
        output[0].append("pump_on")
        output[1].append("bool_")
        output[0].append("btl_fire")
        output[1].append("bool_")
        output[0].append("scan_datetime")
        output[1].append("float64")

        return output

    def _location_fix(self, b1, b2, b3, b4, b5, b6, b7):
        '''
        Function: _location_fix
        Description: This function takes the location coordinates provided by
        the instrument and performs a fix to output in the correct format.
        Parameter(s): self, b1, b2, b3, b4, b5, b6, b7
        Return(s): output
        '''

        lat = (int(b1, 16) * 65536 + int(b2, 16) * 256 + int(b3, 16)) / 50000
        lon = (int(b4, 16) * 65536 + int(b5, 16) * 256 + int(b6, 16)) / 50000

        flag_new_fix = False

        mask_lat_pos = 0x80
        mask_lon_pos = 0x40
        mask_new_fix = 0x01

        if int(b7, 16) & mask_lat_pos:
            lat = lat * -1
        if int(b7, 16) & mask_lon_pos:
            lon = lon * -1
        if int(b7, 16) & mask_new_fix:
            flag_new_fix = True

        output = "{0}, {1}, {2}".format(lat, lon, flag_new_fix)

        return output

    def _reverse_bytes(self, hex_time):
        '''
        Function: _reverse_bytes
        Description: This function takes the hex time and simply reverses it
        and adds a hex prefix to produce a reversed hex time.
        Parameter(s): self, hex_time
        Return(s): reverse_hex
        '''
        time = hex_time.decode("utf-8")
        if re.match("0x", time):
            time = time[2:]
        if(len(time) % 2) == 1:
            time = "0" + time

        # Start flipping
        list_1 = []
        reverse_hex = ""
        for y in range(0, len(time), 2):
            list_1.append(time[y:y + 2])
        for y in range(len(list_1) - 1, -1, -1):
            reverse_hex = reverse_hex + list_1[y]

        # Add prefix to make python hex compatible
        reverse_hex = "0x" + reverse_hex

        return reverse_hex


    def _sbe_time(self, hex_time, sbe_type):
        '''
        Function: _sbe_time
        Description: This function takes the relevant information extracted
        from the xml config file and the hex file to get the data collected
        in the hex file produced by the SBE instrument.
        Parameter(s): self
        Return(s): time
        '''

        #Values as stated by SBE in manual
        scan_start = datetime.datetime(1970, 1, 1, 0, 0, 0)
        nmea_start = datetime.datetime(2000, 1, 1, 0, 0, 0)

        #Catch if 8 characters was not input
        if re.match("0x", hex_time):
            hex_time = hex_time[2:]
        if len(hex_time) > 8:
            raise Exception("Hex string too long to be SBE formatted time: " +
                            hex_time)
        elif len(hex_time) < 8:
            raise Exception("Hex string too short to be SBE formatted time: " +
                            hex_time)

        seconds = datetime.timedelta(seconds=int(hex_time, 16))

        if sbe_type == "scan":
            time = scan_start + seconds
            return time.replace(tzinfo=timezone("UTC")).timestamp()
        elif sbe_type == "nmea":
            time = nmea_start + seconds
            return time.replace(tzinfo=timezone("UTC")).timestamp()
        else:
            raise Exception("Please choose 'nmea' or 'scan' for second input\
            to  _sbe_time()")


    def _sbe_time_create(self, utc_time, sbe_type="scan"):
        '''
        Function: _sbe_time_create
        Description: This function creates a date & time format of when the
        scans took place.
        Parameter(s): self, utc_time, sbe_type
        Return(s): output
        '''
        scan_start = datetime.datetime(1970, 1, 1, 0, 0, 0,
                                       tzinfo=timezone("UTC"))
        nmea_start = datetime.datetime(2000, 1, 1, 0, 0, 0,
                                       tzinfo=timezone("UTC"))
        utc_time = utc_time.replace(tzinfo=timezone("UTC"))
        if sbe_type == "scan":
            time = utc_time - scan_start
        elif sbe_type == "nmea":
            time = utc_time - nmea_start

        hex_time = hex(int(time.total_seconds()))

        output = self._reverse_bytes(bytearray(hex_time, "utf-8"))

        return output

    def _sbe_time_seq(self, data_list):
        '''
        Function: _sbe_time_seq
        Description: This function takes the date time information generated by
        the instrument and formats it.
        Parameter(s): self, data_list
        Return(s): output
        '''

        #Pull out the second to last line of the commands in .hex, then pull
        #pull out the datetime info at the end of the line, then format
        start_month = self.raw_comments[-2][-4]
        start_day = self.raw_comments[-2][-3]
        start_year = self.raw_comments[-2][-2]
        start_time = self.raw_comments[-2][-1]

        start_scan_time = datetime.datetime.strptime(start_month + start_day +
                                                     start_year + start_time,
                                                     "%b%d%Y%H:%M:%S")
        current_scan_time = start_scan_time
        hz_counter = 0
        output = []
        for l in data_list:
            if hz_counter >= 24:
                hz_counter = 0
                current_scan_time = current_scan_time +\
                datetime.timedelta(seconds=1)
            l = l + str(
                current_scan_time.replace(tzinfo=timezone("UTC")).timestamp())
            output.append(l)
            hz_counter += 1
            return output

    def _flag_status(self, flag_char, scan_number):
        '''
        Function: _flag_status
        Description: This function uses the flag value to and the appropriate
        masking to set output to determine a set of conditions of the
        instrument during its data collection process.
        Parameter(s): self, flag_char, scan_number
        Return(s): output
        '''
        mask_pump = 0x1
        mask_bottom_contact = 0x2
        mask_water_sampler_interface = 0x4
        mask_modem_carrier_detect = 0x8

        output = ""

        if flag_char & mask_pump:
            output += "pump on, "
        else:
            output += "pump off, "

        if flag_char & mask_bottom_contact:
            output += "no bottom contact, "
        else:
            output += "bottom contact, "

        if flag_char * mask_water_sampler_interface:
            output += "G.O. 1015 water sampler confirm signal OR manual\
                pump installed, "
        else:
            output += "no signal from G.O. 1015 water sampler OR no manual\
                pump installed, "

        if flag_char & mask_modem_carrier_detect:
            output += "modem carrier signal detected"
        else:
            output += "no modem carrier signal detected"
        irint(line. ordering)

        return output

    def _pump_status(self, flag_char):
        '''
        Function: _pump_status
        Description: This function checks the reported pump status of the
        instrument during data collection and returns true if pump was on, and
        false if it wasn't.
        Parameter(s): self, flag_char
        Return(s): True or False
        '''
        mask_pump = 1

        if int(flag_char) & mask_pump:
            return True
        else:
            return False

    def _bottle_fire(self, flag_char):
        '''
        Function: _bottle_fire
        Description: This function checks whether the instrument fired a bottle
        during data collection and returns true if fired or false if not.
        Parameter(s): self, flag_char
        Return(s): True or False
        '''
        mask_bottle_fire = 0x4

        if int(flag_char, 16) & mask_bottle_fire:
            return True
        else:
            return False

    def _digiquartz_temp_correction(self, hex_in):
        '''
        Function: _digiquartz_temp_correction
        Description: This function corrects the hex temperature into a usable
        format.
        Parameter(s): self, hex_in
        Return(s): T_d
        '''
        charm = format(hex_in, "d")
        T_d = (AD590M * charm) + AD590B
        return T_d

    def _load_hex(self):
        '''
        Function: _load_hex
        Description: This function takes the hex file and processes it for
        parsing by removing irrelevant information and separating comments from
        the raw data.
        Parameter(s): self
        Return(s): none
        '''
        split_lines = self.raw_hex.splitlines()
        self.raw_bytes = [l.strip().encode("utf-8") for l in split_lines if
                          not l.startswith("*")]
        self.raw_comments = [l.strip().split() for l in split_lines if
                             l.startswith("*")]

    def _check_scan_lengths(self):
        '''
        Function: _check_scan_lengths
        Description: This function compares the length of the data in the hex
        file and the scan length as reported by the config to make sure that
        the two files are the right pair.
        Parameter(s): self
        Return(s): none
        '''
        if not all([len(scan) == self.scan_length for scan in
                    self.raw_bytes]):
            raise ValueError("The data length does not match the expected\
            length from the config")

    def calculate_scan_byte_length(self):
        '''
        Function: _calculate_scan_byte_length
        Description: This function calculates the scan byte length to ensure
        that it is the right size before parsing.
        Parameter(s): self
        Return(s): int
        '''
        #Length is extra info that's always present and takes 6 chars
        length = 6
        length += 6 * (5 - self.config["FrequencyChannelsSuppressed"])
        length += 3 * (8 - self.config["VoltageWordsSuppressed"])
        length += 8 * self.config["ScanTimeAdded"]
        length += 8 * self.config["NmeaTimeAdded"]
        length += 6 * self.config["NmeaDepthDataAdded"]
        length += 14 * self.config["NmeaPositionDataAdded"]
        length += 6 * self.config["SurfaceParVoltageAdded"]

        return int(length)

    @classmethod
    def from_paths(cls, raw_hex_path, xml_config_path, encoding="cp437"):
        '''
        Function: from_paths
        Description: This function creates an instance of the SBEReader class
        as an object with the data from the file provided as strings after
        enforcing an encoding for character interpretation.
        Parameter(s): cls, raw_hex_path, xml_config_path, encoding
        Return(s): instance
        '''
        with open(raw_hex_path, encoding=encoding) as raw_hex_file, open(
            xml_config_path, encoding=encoding) as xml_config_file:
            return cls(raw_hex_file.read(), xml_config_file.read())

    @classmethod
    def from_dict(cls, data):
        '''
        Function: _from_dict
        Description: This function creates an instance of the SBE reader class
        as a dict.
        Parameter(s): cls, data
        Return(s): instance
        '''
        instance = cls(data["raw_hex"], data["xml_config"])
        try:
            instance._parsed_scans = data["_parsed_scans"]
        except KeyError:
            pass
        return instance

    @property
    def parsed_scans(self):
        '''
        Function: _parsed_scans
        Description: This function is a getter attribute for the parsing the
        scans from the instrument.
        Parameter(s): self
        Return(s): _parse_scans
        '''
        try:
            return self._parse_scans
        except AttributeError:
            self._parse_scans = pd.concat((self._parse_scans(),
                                           self._parse_scans_meta().reshape(
                                               self._parse_scans_meta().size,
                                               1)), axis=1)
            return self._parse_scans

    @property
    def parsed_config(self):
        '''
        Function: _parsed_config
        Description: This function is a getter attribute that returns the
        parsed config.
        Parameter(s): self
        Return(s): config
        '''
        return self.config

    @property
    def scan_length(self):
        '''
        Function: _scan_length
        Description: This function is a getter attribute for calculating scan
        byte length.
        Parameter(s): self
        Return(s): calculate_scan_byte_length()
        '''
        return self.calculate_scan_byte_length()

    @property
    def frequency_channels(self):
        '''
        Function: frequency_channels
        Description: This function is a getter attribute for
        frequency_channels.
        Parameter(s): self
        Return(s): none
        '''
        pass
