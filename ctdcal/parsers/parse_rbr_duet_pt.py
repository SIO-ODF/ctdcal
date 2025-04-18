"""
Parses RBR duet pressure/temperature raw files in sqlite3 format. Continuous data
are separated and trimmed into casts using known start and end times from ctd data.
"""
import logging
import sqlite3
from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path

import pandas as pd

from ctdcal.common import validate_file

log = logging.getLogger(__name__)


class Parser(object):
    def __init__(self, cast_list, indir, cnvdir, btldir):
        self.indir = indir
        self.cnvdir = cnvdir
        self.btldir = btldir
        self.raw_files = sorted(Path(self.indir).glob('*.rsk'))
        self.raw_file_timestamps = [dt.strptime(ts.stem[-13:], '%Y%m%d_%H%M') for ts in self.raw_files]
        self.ctd_files = [Path(self.cnvdir, '%s.pkl' % cast) for cast in cast_list]
        self.start_times = self.get_start_times()

    def get_start_times(self):
        cast_start_times = dict()
        for ctd_file in self.ctd_files:
            cast_time = pd.read_pickle(ctd_file).iloc[0]['scan_datetime']
            cast_start_times[ctd_file.stem] = dt.fromtimestamp(cast_time)
        return cast_start_times


class Cast(object):
    """
    Container for cast data and metadata with methods identifying the source raw
    data file, and parsing the raw data in sqlite3 format.
    """
    def __init__(self, cast_id):
        self.cast_id = cast_id
        self.raw_file = None
        self.bottle_times = None
        self.data = None

    def load_bottle_times(self, indir, timecol):
        """
        Gets cast bottle fire times from the ctd bottle files. Sets the bottle_times
        attribute as a 1d array or pandas Series.

        Parameters
        ----------
        indir : str or Path-like
            bottle file directory
        timecol : str
            name of the time column
        """
        infile = validate_file(Path(indir, '%s_btl_mean.pkl' % self.cast_id))
        self.bottle_times = pd.read_pickle(infile)[timecol]

    def select_raw_file(self, parser):
        """
        Locates the RBR raw data file which contains the cast based on cast
        start times and the raw file names/timestamps.

        Parameters
        ----------
        parser : Parser object

        Returns
        -------

        """
        # find start time for the cast
        start_time = parser.start_times[self.cast_id]
        # sort through raw files by timestamp range and stop at the one which
        # contains the cast time
        ts = None
        for ts in sorted(parser.raw_file_timestamps, reverse=True):
            if start_time > ts:
                break
        self.raw_file = parser.raw_files[parser.raw_file_timestamps.index(ts) + 1]

    def parse_raw(self, start_time):
        """
        Queries RBR raw data files in sqlite3 format for upcast and pressure calibration
        data. Sets the data attribute to a DataFrame of raw upcast data and corrected
        pressure data.
        TODO: needs some data validation including what happens when the cal interval isn't
            found or whether the calibrated pressure is better than the factory offset-
            corrected pressure.

        Parameters
        ----------
        start_time : datetime.datetime
            upcast start time

        """
        # construct the sql queries
        upcast_start = self.bottle_times.min()
        upcast_stop = self.bottle_times.max() + 300   # add 5 min
        cast_query = 'SELECT * FROM data WHERE tstamp BETWEEN ? AND ?'
        cast_query_params = (upcast_start * 1000, upcast_stop * 1000)
        pcal_start = (start_time - td(minutes=30)).timestamp()
        pcal_stop = (start_time - td(minutes=20)).timestamp()
        cal_query = 'SELECT * FROM data WHERE tstamp BETWEEN ? AND ?'
        cal_query_params = (pcal_start * 1000, pcal_stop * 1000)
        offset_meta_query = 'SELECT value FROM parameterKeys WHERE key="PRESSURE"'

        # execute the queries and create the dataframes
        with sqlite3.connect(self.raw_file) as con:
            data = pd.read_sql_query(cast_query, con, params=cast_query_params)
            cal_pressures = pd.read_sql_query(cal_query, con, params=cal_query_params)['channel02']

            cur = con.cursor().execute(offset_meta_query)
            offset_meta = float(cur.fetchone()[0])
        # convert timestamps
        data['tstamp'] = pd.to_datetime(data['tstamp'], unit='ms')
        # calibrate pressure using pre-cast data (pressure_cor) and using
        # factory offset (sea_pressure)
        data['pressure_cor'] = data['channel02'] - cal_pressures.mean()
        data['sea_pressure'] = data['channel02'] - offset_meta
        # rename columns and set data attribute
        self.data = data.rename(
                columns={'channel01': 'temp_raw', 'channel02': 'pressure_raw', 'channel03': 'ptemp_raw'}
        )


def parse_cast(cast_id, parser):
    """
    Parses a raw RBR sensor data file into a single cast data file which can
    be processed into a reference bottle temperature file.

    Parameters
    ----------
    cast_id : str
        cast name.
    parser : Parser object
        metadata required for the cast parser.

    Returns
    -------
    Cast object

    """
    # initialize the cast
    cast = Cast(cast_id)
    cast.load_bottle_times(parser.btldir, 'scan_datetime')
    # find the source raw file
    try:
        cast.select_raw_file(parser)
    except IndexError:
        log.warning("No RBR duet PT data found for cast %s. Moving along." % cast_id)
        return None
    # parse it
    cast.parse_raw(parser.start_times[cast_id])
    return cast
