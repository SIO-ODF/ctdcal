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

    Input:
    line - a sequence specified by format in _parse_scans_meta

    Output:
    __untitled - a string of converted data

    Format:
    Lat, Lon, bottle fire status, NMEA time, Scan time
    '''
    def _breakdown(line):
        # skip line[0]
        #convert to lat/lon
        _location_fix(line[1], line[2], line[3], line[4], line[5], line[6], line[7])

        '''
        Depth is not implemented yet, SBE does not know how they did it.
        '''
        #skip line[8]

        '''
        NMEA time is LSB, in seconds since January 1, 2000. check it works correctly or break down bytes again.
        '''
        #line[9]
        _sbe_time(_reverse_bytes(line[9]), 'nmea')

        '''
        Pressure temp, status and modulo required by SBE format.
        '''
        #skip line[n+1]

        '''
        Scan time is LSB, 8 char long, in seconds since January 1, 1970
        '''
        #line[13]
        _sbe_time(_reverse_bytes(line[13]), 'scan')

        return None

    def _location_fix(b1, b2, b3, b4, b5, b6, b7):
        """Determine location from SBE format.

        Input:
        b1, b2, b3: three components of Latitude
        b4, b5, b6: three components of Longitude
        b7: sign for lat/lon, is it a new fix

        Output:
        Tuple with three components in order:
        latitude, longitude, new fix?
        Latitude is a float
        Longitude is a float
        If it is a new fix it will be true, otherwise false
        """
        lat = (b1 * 65536 + b2 * 256 + b3)/50000
        lon = (b4 * 65536 + b5 * 256 + b6)/50000

        '''
        If bit 1 in byte_pos is 1, this is a new position
        If bit 8 in byte_pos is 1, lat is negative
        If bit 7 in byte_pos is 1, lon is negative
        '''
        flag_new_fix = False

        mask_lat_pos = 0x80
        mask_lon_pos = 0x40
        mask_new_fix = 0x01

        if b7 & mask_lat_pos:
            lat = lat * -1
        if b7 & mask_lon_pos:
            lon = lon * -1
        if b7 & mask_new_fix:
            flag_new_fix = True

        #print('Latitude: ' + lat)
        #print('Longitude: ' + lon)
        return (lat, lon, flag_new_fix)

    def _reverse_bytes(hex_time):
        """Reverse hex time according to SBE docs.
        Split number by every two chars, then recombine them in reverse order.

        This makes no assumptions about length or whatever. If there are an odd number of characters
        it adds a 0 before the string to pad it out.

        The final output is guaranteed to have '0x' at the beginning.
        """
        time = hex_time
        if re.match('0x', time):
            time = time[2:]
        print(time)
        if (len(time) % 2) == 1:
            time = '0' + time

        """Start flipping."""
        list_1 = []
        reverse_hex = ''
        for y in range(0, len(time), 2):
            list_1.append(time[y:y + 2])
        for y in range(len(list_1) - 1, -1, -1):
            reverse_hex = reverse_hex + list_1[y]
        """Add prefix to make python hex compatible"""
        reverse_hex = '0x' + reverse_hex
        return reverse_hex

    def _sbe_time(hex_time, sbe_type):
        """Convert raw data from sbe .hex file to appropriate datetime, then return string.
        Assumes hex_string has already been reversed by every two chars.
        Needs to be checked for timezone problems in the future.

        Input:
        hex_time: raw hex corresponding to nmea time or scan time. Must be 8 characters.
        sbe_type: either "nmea" or "scan"

        Output:
        A string to be ready printed to csv.

        """
        #values as state by SBE in manual
        scan_start = datetime.datetime(1970, 1, 1, 0, 0, 0)
        nmea_start = datetime.datetime(2000, 1, 1, 0, 0, 0)

        #catch if 8 characters was not input
        if re.match('0x', hex_time):
            hex_time = hex_time[2:]
        if len(hex_time) > 8:
            raise Exception('Hex string too long to be SBE formatted time')
        if len(hex_time) < 8:
            raise Exception('Hex string too short to be SBE formatted time')

        seconds = datetime.timedelta(seconds = int(hex_time, 16))

        if sbe_type == "scan":
            time = scan_start + seconds
            return str(time)
        elif sbe_type == "nmea":
            time = nmea_start + seconds
            return str(time)
        else:
            raise Exception('Please choose "nmea" or "scan" for second input to _sbe_time()')

    def _flag_status(flag_char, scan_number):
        """An attempt to reverse engineer what the single status character means.

        bit 1 = pump status, 1 = on, 0 = off
        bit 2 = bottom contact status, no bottom contact sensor to test with
        bit 4 = bottle fire status, 1 = fired, 0 = not fired
        bit 8 = not used?, default to 0

        bit 2 seems to default to 1 when no bottom contact sensor is installed
        bit 4 when set notes a bottle has been fired, no positional info

        NOT FINISHED
        """

        mask_pump = 0x1
        mask_bottom_contact = 0x2
        mask_bottle_fire = 0x4

        if flag_char & mask_bottle_fire:
            print("Bottle fire at scan: " + str(scan_number))
        if flag_char & mask_pump:
            print("Pump on")

    def _bottle_fire(flag_char):
        """Determine if a scan is around a bottle firing as marked by SBE.
        A simplification of _flag_status.

        bit 4 = bottle fire status, 1 = fired, 0 = not fired
        bit 4 when set notes a bottle has been fired, no positional info

        Input:
        flag_char - a single hex char that

        Output:
        Returns True if scan is a part of a bottle fire,
        or False if scan is not a part of a bottle fire.

        Need to check for off-by-one errors
        """

        mask_bottle_fire = 0x4

        if flag_char & mask_bottle_fire:
            return(True)
        else:
            return(False)

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
