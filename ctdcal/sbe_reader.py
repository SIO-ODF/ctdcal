#!/bin/bash/env
"""
This module is used to parse the xml config and hex file provided by user. It
takes the xml config and extracts the useful settings used during a CTD cast,
and uses that to check and parse the hex file.
"""

import datetime
import re
import struct
import xml.etree.cElementTree as ET

import numpy as np
from pytz import timezone


class SBEReader:
    """Code originally written by Andrew Barna, January-March 2016."""

    def __init__(self, raw_hex, xml_config):
        """expects long character string inputs"""
        self.raw_hex = raw_hex
        self.xml_config = xml_config
        self._parse_config()
        self._load_hex()
        self._check_scan_lengths()

    def _load_hex(self):
        split_lines = self.raw_hex.splitlines()
        self.raw_bytes = [
            line.strip().encode("utf-8")
            for line in split_lines
            if not line.startswith("*")
        ]
        # next few lines are to grab start_scan_time
        self.raw_comments = [
            line.strip().split() for line in split_lines if line.startswith("*")
        ]

    def _check_scan_lengths(self):
        if not all([len(scan) == self.scan_length for scan in self.raw_bytes]):
            raise ValueError(
                "The data length does not match the expected length from the config"
            )

    def _parse_scans(self):
        """The order according to the SBE docs are:
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
        """

        # parse the frequencies
        # scan_list = list(scan)
        num_frequencies = 5 - self.config["FrequencyChannelsSuppressed"]
        num_voltages = 8 - self.config["VoltageWordsSuppressed"]
        flag_spar = 0

        if self.config["SurfaceParVoltageAdded"]:
            flag_spar = 1

        the_bytes = b"".join(self.raw_bytes)

        # specify how each line of raw data should be broken down
        unpack_str = (
            "6s" * num_frequencies
            + "3s" * num_voltages
            + "2s" * flag_spar
            + "4s" * flag_spar
            + "{}s".format(
                self.scan_length
                - num_voltages * 3
                - num_frequencies * 6
                - flag_spar * 6
            )
        )
        measurements = np.array(
            [
                [int(x, 16) for x in line]
                for line in struct.iter_unpack(unpack_str, the_bytes)
            ]
        )

        # this uses numpy magic to just "apply" the needed operations to the correct part of the array all in one go
        measurements[:, :num_frequencies] = measurements[:, :num_frequencies] / 256
        measurements[:, num_frequencies : num_frequencies + num_voltages] = 5 * (
            1
            - (measurements[:, num_frequencies : num_frequencies + num_voltages] / 4095)
        )
        measurements = measurements[
            :, 0:-1
        ]  # remove after applying _parse_scans_meta in here
        # print(measurements.shape)
        return measurements

    def _parse_scans_meta(self):
        """The order according to the SBE docs are:
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
        """

        # parse the frequencies
        # scan_list = list(scan)
        num_frequencies = 5 - self.config["FrequencyChannelsSuppressed"]
        num_voltages = 8 - self.config["VoltageWordsSuppressed"]

        string_order = ["scan"]
        # flags to determine how many bytes to break by
        flag_nmea_pos = 0
        flag_nmea_depth = 0
        flag_nmea_time = 0
        flag_scan_time = 0
        flag_spar = 0
        # put in to make code readable later on
        flag_pressure_temp = 1
        flag_ctd_status = 1
        flag_modulo = 1
        if self.config["SurfaceParVoltageAdded"]:
            flag_spar = 1
            # string_order.append('spar') #already added into mass line, does not need additional flag
        if self.config["NmeaPositionDataAdded"]:
            flag_nmea_pos = 1
            string_order.append("nmea_pos")
            # Depth is here for completeness but not implemented,
            # after email chain showed SBE no longer knows how they did it.
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

        # specify how each line of raw data should be broken down
        # this format is fed directly into _breakdown, so it needs to be tracked there
        unpack_str = (
            str(num_frequencies * 6 + num_voltages * 3 + flag_spar * 6)
            + "s"
            + "14s" * flag_nmea_pos
            + "6s" * flag_nmea_depth
            + "8s" * flag_nmea_time
            + "3s" * flag_pressure_temp
            + "1s" * flag_ctd_status
            + "2s" * flag_modulo
            + "8s" * flag_scan_time
        )
        # breaks down line according to columns specified by unpack_str and converts from hex to decimal MSB first
        # needs to be adjusted for LSB fields (NMEA time, scan time), break into two lines? move int(x,16) outside
        measurements = [line for line in struct.iter_unpack(unpack_str, the_bytes)]
        # print(unpack_str)
        # measurements_2 = np.array([self._breakdown(line, string_order) for line in measurements])
        measurements_2 = [self._breakdown(line, string_order) for line in measurements]
        # if no time is enabled, fake the scan timestamp from info in the .hex file
        if flag_scan_time == 0:
            measurements_3 = self._sbe_time_seq(measurements_2)
            #            print(True)
            measurements_3 = np.array(measurements_3)
            return measurements_3
        # print(measurements_3)
        return np.array(measurements_2)

    def _breakdown(self, line, ordering):
        """
        Convert hex metadata to science useable data.

        Steps:
        1. Throw away away first part of the line
        2. Start processing the rest of the line
        3. Return a sequence (tuple? list?)

        Input:
        line - a sequence specified by format in _parse_scans_meta
        ordering - the order the string comes in

        Output:
        output - a string of converted data according to format

        Format:
        Lat, Lon, pressure temp, bottle fire status, NMEA time, Scan time
        """
        output = ""
        # print(ordering, line)
        for x, y in zip(ordering, line):
            if x == "nmea_pos":
                tokens = []
                # print(y)
                for t in struct.iter_unpack("2s2s2s2s2s2s2s", y):
                    for tt in t:
                        tokens.append(tt)
                output = (
                    output
                    + str(
                        self._location_fix(
                            tokens[0],
                            tokens[1],
                            tokens[2],
                            tokens[3],
                            tokens[4],
                            tokens[5],
                            tokens[6],
                        )
                    )
                    + ","
                )
            elif x == "nmea_time":
                output = (
                    output + str(self._sbe_time(self._reverse_bytes(y), "nmea")) + ","
                )
            elif x == "scan_time":
                output = output + str(self._sbe_time(self._reverse_bytes(y), "scan"))
            elif x == "flag_ctd_status":
                output = (
                    output
                    + str(self._pump_status(y))
                    + ","
                    + str(self._bottle_fire(y))
                    + ","
                )  # need to make final comma dynamic later so less errors down the line
            elif x == "pressure_temp":
                output = output + str(int(y, 16)) + ","

        return output

    def _breakdown_header(self):
        """Creates header for metadata. Arrays below are what is expected.

        ['GPSLAT', 'GPSLON', 'new_fix', 'nmea_datetime', 'pressure_temp_int', 'btl_fire', 'scan_datetime'],
        ['float64', 'float64', 'bool_', 'datetime64', 'int_', 'bool_', 'datetime64']
        """
        # temp fix, need to adjust to take in file to adjust as wanted?
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
        # if self.config["ScanTimeAdded"]: #need to fix flagging part
        # output[0].append('scan_datetime')
        # output[1].append('float64')
        output[0].append("scan_datetime")
        output[1].append("float64")

        return output

    def _location_fix(self, b1, b2, b3, b4, b5, b6, b7):
        """Determine location from SBE format.

        Input:
        b1, b2, b3: three components of Latitude
        b4, b5, b6: three components of Longitude
        b7: sign for lat/lon, is it a new fix

        Output:
        Tuple with three components in order:
        latitude, longitude, new fix
        Latitude is a float
        Longitude is a float
        If it is a new fix it will be true, otherwise false
        """
        lat = (int(b1, 16) * 65536 + int(b2, 16) * 256 + int(b3, 16)) / 50000
        lon = (int(b4, 16) * 65536 + int(b5, 16) * 256 + int(b6, 16)) / 50000

        """
        If bit 1 in byte_pos is 1, this is a new position
        If bit 8 in byte_pos is 1, lat is negative
        If bit 7 in byte_pos is 1, lon is negative
        """
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

        output = "{0},{1},{2}".format(lat, lon, flag_new_fix)
        return output

    def _reverse_bytes(self, hex_time):
        """Reverse hex time according to SBE docs.
        Split number by every two chars, then recombine them in reverse order.

        This makes no assumptions about length or whatever. If there are an odd number of characters
        it adds a 0 before the string to pad it out.

        The final output is guaranteed to have '0x' at the beginning.
        """
        time = hex_time.decode("utf-8")
        if re.match("0x", time):
            time = time[2:]
        if (len(time) % 2) == 1:
            time = "0" + time

        """Start flipping."""
        list_1 = []
        reverse_hex = ""
        for y in range(0, len(time), 2):
            list_1.append(time[y : y + 2])
        for y in range(len(list_1) - 1, -1, -1):
            reverse_hex = reverse_hex + list_1[y]
        """Add prefix to make python hex compatible"""
        reverse_hex = "0x" + reverse_hex
        # print('Output of reverse_hex: ' + reverse_hex)
        return reverse_hex

    def _sbe_time(self, hex_time, sbe_type):
        """Convert raw data from sbe .hex file to appropriate datetime, then return string.
        Assumes hex_string has already been reversed by every two chars.
        Needs to be checked for timezone problems in the future.

        Input:
        hex_time: raw hex corresponding to nmea time or scan time. Must be 8 characters.
        sbe_type: either "nmea" or "scan"

        Output:
        A string to be ready printed to csv.

        """
        # values as state by SBE in manual
        scan_start = datetime.datetime(1970, 1, 1, 0, 0, 0)
        nmea_start = datetime.datetime(2000, 1, 1, 0, 0, 0)

        # catch if 8 characters was not input
        if re.match("0x", hex_time):
            hex_time = hex_time[2:]
        if len(hex_time) > 8:
            raise Exception("Hex string too long to be SBE formatted time: " + hex_time)
        if len(hex_time) < 8:
            raise Exception(
                "Hex string too short to be SBE formatted time: " + hex_time
            )

        seconds = datetime.timedelta(seconds=int(hex_time, 16))
        # changed to epoch time at request of analysts
        # epoch time assumes it is UTC all way through
        if sbe_type == "scan":
            # time = scan_start + seconds
            # return time.isoformat()
            time = scan_start + seconds
            return time.replace(tzinfo=timezone("UTC")).timestamp()
        elif sbe_type == "nmea":
            # time = nmea_start + seconds
            # return time.isoformat()
            time = nmea_start + seconds
            return time.replace(tzinfo=timezone("UTC")).timestamp()
        else:
            raise Exception(
                'Please choose "nmea" or "scan" for second input to _sbe_time()'
            )

    def _sbe_time_create(self, utc_time, sbe_type="scan"):
        """Reverse of _sbe_time, create a sbe timestamp in hex to append to file.
        Useful when SBE acquisition software does not have NMEA time or scan time activated.

        The inverse of _sbe_time

        Input:
        utc_time: datetime
        """

        # values as state by SBE in manual
        scan_start = datetime.datetime(1970, 1, 1, 0, 0, 0, tzinfo=timezone("UTC"))
        nmea_start = datetime.datetime(2000, 1, 1, 0, 0, 0, tzinfo=timezone("UTC"))
        utc_time = utc_time.replace(tzinfo=timezone("UTC"))
        if sbe_type == "scan":
            time = utc_time - scan_start
        elif sbe_type == "nmea":
            time = utc_time - nmea_start

        # timedelta method to return seconds from any timedelta format
        hex_time = hex(int(time.total_seconds()))

        output = self._reverse_bytes(bytearray(hex_time, "utf-8"))
        return output

    def _sbe_time_seq(self, data_list):
        """Recreates the scan timestamp if the option was not enabled in SBE acq.
        Accurate to 1 second/24hz, as it uses the start time in the second to last line of the .hex file.

        Sequencer will run over a processed file/list? and append the scan time at the end.
        """
        # Pull out the second to last line of the comments in .hex,
        # then pull out the datetime info at the end of the line, then format
        start_month = self.raw_comments[-2][-4]
        start_day = self.raw_comments[-2][-3]
        start_year = self.raw_comments[-2][-2]
        start_time = self.raw_comments[-2][-1]

        start_scan_time = datetime.datetime.strptime(
            start_month + start_day + start_year + start_time, "%b%d%Y%H:%M:%S"
        )
        current_scan_time = start_scan_time
        hz_counter = 0
        output = []
        for line in data_list:
            if hz_counter >= 24:
                hz_counter = 0
                current_scan_time = current_scan_time + datetime.timedelta(seconds=1)
            line = line + str(
                current_scan_time.replace(tzinfo=timezone("UTC")).timestamp()
            )
            output.append(line)
            hz_counter += 1
        return output

    def _flag_status(self, flag_char, scan_number):
        """Decode SBE flag bit, as referenced on SBE 11pV2, pg 66.
        CTD status:
            Bit 0 Pump status- 1=pump on, 0=pump off.
            Bit 1 Bottom contact switch status -
                1 = switch open (no contact), 0 = switch closed.
            Bit 2 G.O. 1015 water sampler interface confirm signal or
                manual pump control signal -
                1 = Deck Unit detects confirm signal from G.O. 1015 or detects manual pump control installed in 9plus,
                0 = not detected.
            Bit 3 CTD modem carrier detect -
                0 = CTD modem detects Deck Unit modem carrier signal, 1 = not detected.

        Output could be cleaned up if it ever needs to be used.
        """

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
        if flag_char & mask_water_sampler_interface:
            output += (
                "G.O. 1015 water sampler confirm signal OR manual pump installed, "
            )
        else:
            output += "no signal from G.O. 1015 water sampler OR no manual pump, "
        if flag_char & mask_modem_carrier_detect:
            output += "modem carrier signal detected"
        else:
            output += "no modem carrier signal detected"

        return output

    def _pump_status(self, flag_char):
        """Decode SBE flag bit, as referenced on SBE 11pV2, pg 66.
        CTD status:
            Bit 0 Pump status- 1=pump on, 0=pump off.
        """
        mask_pump = 1
        if bytes.decode(flag_char).isnumeric() is False:
            flag_char = b"3"
        if int(flag_char) & mask_pump:
            return True
        else:
            return False

    def _bottle_fire(self, flag_char):
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

        if int(flag_char, 16) & mask_bottle_fire:
            return True
        else:
            return False

        return None

    # no calls to this, what are AD590M and AD590B? calib coefs prob
    # def _digiquartz_temp_correction(self, hex_in):
    #     """Digiquartz pressure sensor temperature correction for internal probe."""
    #     charm = format(hex_in, "d")
    #     T_d = (AD590M * charm) + AD590B
    #     return T_d

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
        sensors = {}

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
                self.config[key] = bool(
                    int(config.find("./Instrument/{}".format(key)).text)
                )
            except AttributeError:
                raise AttributeError("Could not find {} in XMLCONF".format(key))
            except ValueError:
                raise ValueError("{} Value is not truthy".format(key))

        """Pokedex is a dict of {Sensor index numbers from the config:sensor info}
        Assume that sensor index number is the order the sensors have been entered into the file.
        Therefore, it will be Frequency instruments first, then Voltage instruments.
        Call their index number (starting at 0) in order to pull out the info.

        """
        # pokedex = {}
        for x in config.iter("Sensor"):
            # print(ET.tostring(x))
            """Start creating single sensor dictionary."""
            bulbasaur = {}
            bulbasaur["SensorID"] = x.attrib["SensorID"]
            # load all values into dict - beware of overwriting #NEED TO FIX
            for children in x:
                for y in children.iter():
                    try:
                        bulbasaur[y.tag] = float(y.text)
                    except (TypeError, ValueError):
                        bulbasaur[y.tag] = (
                            str(y.text).replace("\n", "").replace(" ", "")
                        )

                """Add sensor to big dictionary."""
                # pokedex[x.attrib['index']] = bulbasaur
                sensors[int(x.attrib["index"])] = bulbasaur
        # sensors.append(pokedex)
        self.config["Sensors"] = sensors

    @classmethod
    def from_paths(cls, raw_hex_path, xml_config_path, encoding="cp437"):
        with open(raw_hex_path, encoding=encoding) as raw_hex_file, open(
            xml_config_path, encoding=encoding
        ) as xml_config_file:
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
        length = 6  # Extra info alway present takes 6 "chars"
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
        """
        Wrapper for _parse_scans and _parse_scans_meta. Returns np.ndarray
        """
        # try:
        #     return self._parse_scans
        # except AttributeError:
        #     self._parse_scans = np.concatenate((self._parse_scans(), self._parse_scans_meta().reshape(self._parse_scans_meta().size,1)), axis = 1)
        #     return self._parse_scans
        return np.concatenate(
            (
                self._parse_scans(),
                self._parse_scans_meta().reshape(self._parse_scans_meta().size, 1),
            ),
            axis=1,
        )
        # return self._parse_scans

    def parsed_config(self):
        return self.config

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
