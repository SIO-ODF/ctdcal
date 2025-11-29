"""
Parses RBR duet pressure/temperature raw files in sqlite3 format. Continuous data
are separated and trimmed into casts using known start and end times from ctd data.
"""
import logging
import sqlite3
from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import correlate

from ctdcal.common import validate_file

log = logging.getLogger(__name__)


class Parser(object):
    def __init__(self, cast_list, indir, cnvdir, btldir):
        self.indir = indir
        self.cnvdir = cnvdir
        self.btldir = btldir
        self.raw_files = sorted(Path(self.indir).glob('*.rsk'))
        self.ctd_files = [Path(self.cnvdir, '%s.pkl' % cast) for cast in cast_list]
        self.raw_file_end_times = {}
        self.meta = {}
        self.start_times = self.get_start_times()
        self.get_raw_file_metadata()

    def get_raw_file_metadata(self):
        # query the db files
        t_range_query = 'SELECT MIN(tstamp), MAX(tstamp) FROM data'
        clock_drift_query = 'SELECT loggerTimeDrift FROM deployments'
        for raw_file in self.raw_files:
            coeffs = dict()
            with sqlite3.connect(raw_file) as con:
                # query for timestamp start and end
                cur = con.cursor().execute(t_range_query)
                coeffs['t_range'] = [float(t) for t in cur.fetchone()]
                # query for clock drift
                cur = con.cursor().execute(clock_drift_query)
                coeffs['t_drift'] = int(cur.fetchone()[0])
            # convert end time to utc datetime
            self.raw_file_end_times[raw_file.stem] = coeffs['t_range'][1] / 1000
            # set coeffs attribute for the raw file
            self.meta[raw_file.stem] = coeffs

    def get_start_times(self):
        cast_start_times = dict()
        for ctd_file in self.ctd_files:
            cast_start_times[ctd_file.stem] = float(pd.read_pickle(ctd_file).iloc[0]['scan_datetime'])
        return cast_start_times


class Cast(object):
    """
    Container for cast data and metadata with methods identifying the source raw
    data file, and parsing the raw data in sqlite3 format.
    """
    def __init__(self, cast_id):
        self.cast_id = cast_id
        self.raw_file = None
        # self.start_time = None
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
        for raw_file in sorted(parser.raw_files):
            if start_time < parser.raw_file_end_times[raw_file.stem]:
                self.raw_file = raw_file
                return

    def parse_raw(self, start_time):
        """
        Queries RBR raw data files in sqlite3 format for upcast and pressure calibration
        data. Sets the data attribute to a DataFrame of raw upcast data and corrected
        pressure data.

        Parameters
        ----------
        start_time : float
            upcast start time

        """
        # TODO: needs some data validation including what happens when the cal interval isn't
        #     found or whether the calibrated pressure is better than the factory offset-
        #     corrected pressure.
        # construct the sql queries
        upcast_start = self.bottle_times.min()
        upcast_stop = self.bottle_times.max() + 300   # add 5 min
        cast_query = 'SELECT * FROM data WHERE tstamp BETWEEN ? AND ?'
        # cast_query_params = (upcast_start * 1000, upcast_stop * 1000)
        # ## selecting the whole cast for troubleshooting 2025-04-22
        cast_query_params = (start_time * 1000, upcast_stop * 1000)
        pcal_start = start_time - 1800
        pcal_stop = start_time - 1200
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

    def synchronize_cast_data(self, cnvdir):
        ctdfile = Path(cnvdir, '%s.pkl' % self.cast_id)
        reference = pd.read_pickle(ctdfile)[['CTDPRS', 'scan_datetime']]
        reference['tstamp'] = pd.to_datetime(reference['scan_datetime'], unit='s')
        reference = reference.set_index('tstamp')['CTDPRS'].resample('1s').mean()

        measured = self.data[['sea_pressure', 'tstamp']]
        # measured['tstamp'] = pd.to_datetime(measured['tstamp'])
        measured = measured.set_index('tstamp')['sea_pressure'].resample('1s').mean()
        #
        t1 = max(measured.index[0], reference.index[0])
        t2 = min(measured.index[-1], reference.index[-1])
        reference = reference[t1:t2]
        measured = measured[t1:t2]

        # below from https://stackoverflow.com/a/56432463
        n = len(measured)
        sr = 1  # sample rate in sec
        corr = correlate(reference, measured, mode='same', method='direct') / np.sqrt(
                correlate(measured, measured, mode='same', method='direct')[int(n / 2)] *
                correlate(reference, reference, mode='same', method='direct')[int(n / 2)]
        )
        delay_arr = np.linspace(-0.5 * n / sr, 0.5 * n / sr, n)
        delay = round(delay_arr[np.argmax(corr)], 2)

        if abs(delay) > 30:
            log.warning(
                    "Derived time lag for cast %s seems excessive. Leaving this one uncorrected." % self.cast_id
            )
        else:
            log.info("Correcting %s seconds of lag for cast %s." % (delay, self.cast_id))
            self.data['tstamp'] += td(seconds=delay)


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
    cast.select_raw_file(parser)
    if cast.raw_file is None:
        log.warning("No RBR duet PT data found for cast %s. Moving along." % cast_id)
        return None
    # parse it
    cast.parse_raw(parser.start_times[cast_id])
    return cast
