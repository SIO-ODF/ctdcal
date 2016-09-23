#%pylab inline
import xml.etree.ElementTree as ET
import struct
import numpy as np
import re
import sys

#not currently used
def parse_frequency_hex(s):
    #return int(s[0:2], 16) * 256 + int(s[2:4], 16) + int(s[4:6], 16)/256
    return int(s, 16) / 256

#Not currently used
def parse_volt_hex(s):
    return (5 * (1 - (int(s, 16) / 4095)))

class SBEReader():
    """Code originally written by Andrew Barna, January-March 2016."""
    def __init__(self, raw_hex, xml_config):
        '''expects long character string inputs'''
        self.raw_hex = raw_hex
        self.xml_config = xml_config
        self._parse_config()
        self._load_hex()
        self._check_scan_lengths()

    def _load_hex(self):
        split_lines = self.raw_hex.splitlines()
        self.raw_bytes = [l.strip().encode("utf-8") for l in split_lines if not l.startswith("*")]

    def _check_scan_lengths(self):
        if not all([len(scan) == self.scan_length for scan in self.raw_bytes]):
            raise ValueError("The data length does not match the expected length from the config")

    def _parse_scans(self):
        '''The order according to the SBE docs are:
        1) Data from the instrument
          a) Frequency (3 bytes each)
          b) Voltage (12 bits each)
        2) Surface Par (3 bytes)
        3) NMEA lat/lon (7 bytes)
        4) NMEA depth (3 bytes)
        5) NMEA time (4 bytes) (low byte first)
        6) Additional Data from the instrument
          a) Pressure temp (12 bits)
          b) pump status (4 bits)
          c) modulo byte (1 byte)
        7) System time (4 bytes) (low byte first)
        If any of the above are omitted, the length of the hex will be smaller.
        '''

        # parse the frequencies
        #scan_list = list(scan)
        num_frequencies = 5 - self.config["FrequencyChannelsSuppressed"]
        num_voltages = 8 - self.config["VoltageWordsSuppressed"]
        flag_spar = 0

        if self.config["SurfaceParVoltageAdded"]:
            flag_spar = 1

        the_bytes = b"".join(self.raw_bytes)

        #specify how each line of raw data should be broken down
        unpack_str = "6s" * num_frequencies + "3s" * num_voltages + "2s" * flag_spar + "4s" * flag_spar + "{}s".format(self.scan_length - num_voltages * 3 - num_frequencies * 6 - flag_spar * 6)
        measurements = np.array([[int(x, 16) for x in line] for line in struct.iter_unpack(unpack_str, the_bytes)])

        #this uses numpy magic to just "apply" the needed operations to the correct part of the array all in one go
        measurements[:,:num_frequencies] = measurements[:,:num_frequencies] / 256
        measurements[:,num_frequencies:num_frequencies+num_voltages] = (5 * (1 - (measurements[:,num_frequencies:num_frequencies+num_voltages] / 4095)))
        return measurements

    def _parse_scans_meta(self):
        '''The order according to the SBE docs are:
        1) Data from the instrument
          a) Frequency (3 bytes each)
          b) Voltage (12 bits each)
        2) Surface Par (3 bytes)
        3) NMEA lat/lon (7 bytes)
        4) NMEA depth (3 bytes)
        5) NMEA time (4 bytes) (low byte first)
        6) Additional Data from the instrument
          a) Pressure temp (12 bits)
          b) pump status (4 bits)
          c) modulo byte (1 byte)
        7) System time (4 bytes) (low byte first)
        If any of the above are omitted, the length of the hex will be smaller.
        '''

        # parse the frequencies
        #scan_list = list(scan)
        num_frequencies = 5 - self.config["FrequencyChannelsSuppressed"]
        num_voltages = 8 - self.config["VoltageWordsSuppressed"]

        #flags to determine how many bytes to break by
        flag_nmea_pos = 0
        flag_nmea_depth = 0
        flag_nmea_time = 0
        flag_scan_time = 0
        flag_spar = 0
        #put in to make code readable later on
        flag_pressure_temp = 1
        flag_status = 1
        flag_modulo = 1
        if self.config["NmeaPositionDataAdded"]:
            flag_nmea_pos = 1
        if self.config["NmeaDepthDataAdded"]:
            flag_nmea_depth = 1
        if self.config["NmeaTimeAdded"]:
            flag_nmea_time = 1
        if self.config["ScanTimeAdded"]:
            flag_scan_time = 1
        if self.config["SurfaceParVoltageAdded"]:
            flag_spar = 1

        the_bytes = b"".join(self.raw_bytes)

        #specify how each line of raw data should be broken down
        #this format is fed directly into _breakdown, so it needs to be tracked there
        unpack_str = (str(num_frequencies * 6 + num_voltages * 3 + flag_spar * 6) + "s" +
            "2s2s2s" * flag_nmea_pos + "2s2s2s" * flag_nmea_pos + "2s" + flag_nmea_pos +
            "6s" * flag_nmea_depth +
            "8s" * flag_nmea_time +
            "3s" * flag_pressure_temp + "1s" * flag_status + "2s" * flag_modulo + #may be excessive
            "8s" * flag_scan_time)
        #breaks down line according to columns specified by unpack_str and converts from hex to decimal MSB first
        #needs to be adjusted for LSB fields (NMEA time, scan time), break into two lines? move int(x,16) outside
        measurements = [line for line in struct.iter_unpack(unpack_str, the_bytes)]
        measurements_2 = [_breakdown(line) for line in measurements]

        return measurements_2

    '''
    Breakdown the NMEA strings and internal scan time. SBE is shit at this.

    Assumes a tuple is being pushed in, length 14.

    Steps:
    1. Throw away away first part of the line
    2. Start processing the rest of the line
    3. Return a sequence (tuple? list?)

    Reduce magic numbers? Possible?
    '''
    def _breakdown(line):
        # skip line[0]
        #convert to lat/lon
        lat = (line[1] * 65536 + line[2] * 256 + line[3])/50000
        lon = (line[4] * 65536 + line[5] * 256 + line[6])/50000

        '''
        Are these flags even needed? Or just the new fix and masks?
        ----
        If bit 1 in byte_pos is 1, this is a new position
        If bit 8 in byte_pos is 1, lat is negative
        If bit 7 in byte_pos is 1, lon is negative

        Code is written such that all values are 0 at start, and 1 after being flipped.
        '''
        flag_lat_neg = 0
        flag_lon_neg = 0
        flag_new_fix = 0

        mask_lat_pos = 0x80
        mask_lon_pos = 0x40
        mask_new_fix = 0x01

        if line[7] & mask_lat_pos:
            flag_lat_neg = 1
            lat = lat * -1
        if line[7] & mask_lon_pos:
            flag_lon_neg = 1
            lon = lon * -1
        if line[7] & mask_new_fix:
            flag_new_fix = 1

        print('Latitude: ' + lat)
        print('Longitude: ' + lon)

        '''
        Depth is not implemented yet, SBE does not know how they did it.
        '''
        #skip line[8]

        '''
        NMEA time is LSB, in seconds since January 1, 2000. check it works correctly or break down bytes again.
        '''
        #fill in here

        '''
        Pressure temp, status and modulo required by SBE format.
        '''
        #skip line[n+1]

        '''
        Scan time is LSB, 8 char long, in seconds since January 1, 1970
        '''
        #fill in here

        return None

    """Only for values in bools and numeric. """
    def _parse_config(self):
        bools = [
            "SurfaceParVoltageAdded",
            "NmeaPositionDataAdded",
            "NmeaDepthDataAdded",
            "NmeaTimeAdded",
            "ScanTimeAdded",
        ]
        numeric = [
            "FrequencyChannelsSuppressed",
            "VoltageWordsSuppressed",
            ]

        config = ET.fromstring(self.xml_config)
        self.config = {}

        for key in numeric:
            try:
                self.config[key] = int(config.find("./Instrument/{}".format(key)).text)
            except AttributeError:
                raise AttributeError("Could not find {} in XMLCONF".format(key))
            except ValueError:
                raise ValueError("{} Value is not a number".format(key))

        for key in bools:
            try:
                self.config[key] = bool(int(config.find("./Instrument/{}".format(key)).text))
            except AttributeError:
                raise AttributeError("Could not find {} in XMLCONF".format(key))
            except ValueError:
                raise ValueError("{} Value is not truthy".format(key))

    @classmethod
    def from_paths(cls, raw_hex_path, xml_config_path, encoding="cp437"):
        with open(raw_hex_path, encoding=encoding) as raw_hex_file, open(xml_config_path, encoding=encoding) as xml_config_file:
            return cls(raw_hex_file.read(), xml_config_file.read())

    @property
    def scan_length(self):
        return self.calculate_scan_byte_length()

    @property
    def frequency_channels(self):
        pass

    def voltage_channels(self):
        pass

    def calculate_scan_byte_length(self):
        length = 6 # Extra info alway present takes 6 "chars"
        length += 6 * (5 - self.config["FrequencyChannelsSuppressed"])
        length += 3 * (8 - self.config["VoltageWordsSuppressed"])
        length += 8 * self.config["ScanTimeAdded"]
        length += 8 * self.config["NmeaTimeAdded"]
        length += 6 * self.config["NmeaDepthDataAdded"]
        length += 14 * self.config["NmeaPositionDataAdded"]
        length += 6 * self.config["SurfaceParVoltageAdded"]
        return int(length)

    @property
    def parsed_scans(self):
        try:
            return self._parsed_scans
        except AttributeError:
            self._parsed_scans = self._parse_scans()
            return self._parsed_scans

    def to_dict(self, parse_cache=True):
        return {
            "raw_hex": self.raw_hex,
            "xml_config": self.xml_config,
            "_parsed_scans": self.parsed_scans,
            }

    @classmethod
    def from_dict(cls, data):
        instance = cls(data["raw_hex"], data["xml_config"])
        try:
            instance._parsed_scans = data["_parsed_scans"]
        except KeyError:
            pass
        return instance
