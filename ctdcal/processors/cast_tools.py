"""
Tools for cleaning up converted and processed cast data.
"""
import logging
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import signal as sig


log = logging.getLogger(__name__)


class Cast(object):
    """
    Cast data container with methods for separating upcast, filtering
    and trimming soak period.

    Attributes
    ----------
    cast_id : str
        Cast identifier.
    datadir : str or Path-like
        Data directory.
    p_col : str
        Name of pressure column.
    proc : DataFrame
        Pre-processed cast data.
    downcast : DataFrame
        Downast data trimmed from full cast data.
    filtered : array-like
        Filtered cast data.
    trimmed : DataFrame
        Upcast data trimmed of initial soak data.
    ondeck_trimmed : DataFrame
        Full cast data trimmed of on-deck intervals.
    """
    def __init__(self, cast_id, datadir):
        self.cast_id = cast_id
        self.datadir = datadir
        self.p_col = None
        self.proc = None
        self.downcast = None
        self.filtered = None
        self.trimmed = None
        self.ondeck_trimmed = None
        self.load_cast()

    def load_cast(self):
        """
        Read the processed data into a dataframe.
        """
        f = Path(self.datadir, 'converted/%s.pkl' % self.cast_id)
        self.proc = pd.read_pickle(f)

    def parse_downcast(self, data):
        """
        Trim away the upcast starting from the max value of the pressure
        column.

        Parameters
        ----------
        data : DataFrame
            Cast data.
        """
        if self.p_col is None:
            raise AttributeError("Pressure column 'p_col' attribute is not set")
        downcast = data.loc[:self.proc[self.p_col].idxmax()].copy()
        self.downcast = downcast

    def filter(self, data, win_size, win_type="triangle", cols=None):
        """
        Filter processed CTD data using one of three window types (boxcar,
        hann, triangle).

        Previously was 'raw_ctd_filter' in ctdcal.process_ctd.

        Parameters
        ----------
        data : array-like
            Cast data.
        win_size : int
            Window size in number of samples.
        win_type : str
            Window type.
        cols : list
            Column names to filter.
        """
        data = data.copy()
        if win_type == "boxcar":
            win = sig.windows.boxcar(win_size)
        elif win_type == "hann":
            win = sig.windows.hann(win_size)
        elif win_type == "triangle":
            win = sig.windows.triang(win_size)
        else:
            raise AttributeError('Error filtering cast %s! No filter of type %s found.' % (self.cast_id, win_type))
        for col in cols:
            data[col] = sig.convolve(data[col], win, mode="same") / np.sum(win)
        self.filtered = data

    def trim_soak(self, data, win_size, max_soak=20):
        """
        Detects the end of the soak period by identifying local peaks in the
        cast data, and choosing the last one above the threshold.

        data : DataFrame
            Downcast data.
        win_size : int
            Window size in number of samples.
        max_soak : int
            Soak threshold in cast pressure units.
        """
        if self.p_col is None:
            raise AttributeError("Pressure column 'p_col' attribute is not set")
        if self.downcast is None:
            raise AttributeError("'downcast' attribute is not set.")
        if win_size < 0:
            win_size = 0  # must be positive or zero

        # Adapted from: https://stackoverflow.com/questions/48023982/pandas-finding-local-max-and-min
        # Find local peaks
        local_min_vals = data.loc[data[self.p_col] == data[self.p_col].rolling(win_size, center=True).min(),
                                  [self.p_col]]

        # apply max_soak threshold to filter false positives...
        i = local_min_vals.index[-1]
        while local_min_vals.loc[i, self.p_col] > max_soak:
            try:
                local_min_vals.drop(i, inplace=True)
            except KeyError:
                log.warning('Whoa, trouble finding the soak on cast %s! Nothing was trimmed!' % self.cast_id)
                return
            i = local_min_vals.index[-1]

        soak_end = local_min_vals.index[-1]
        self.trimmed = self.downcast.loc[soak_end:]

    def get_details(self):
        """
        Collect cast details for export (see Notes).

        Previously was 'cast_details' in ctdcal.process_ctd.

        Notes
        -----
        The following variables are returned in an array:
        cast_id : cast identifier
        start_time : Time at start of cast (in unix epoch time)
        end_time : Time at end of cast (in unix epoch time)
        bottom_time : Time at bottom of cast (in unix epoch time)
        start_pressure : Pressure at which cast started
        max_pressure : Bottom of the cast pressure
        latitude : Latitude at bottom of cast
        longitude : Longitude at bottom of cast
        altimeter_bottom : Altimeter reading at bottom of cast

        Returns
        -------
        DataFrame
        """
        if self.p_col is None:
            raise AttributeError("Pressure column 'p_col' attribute is not set")
        if self.trimmed is None:
            raise AttributeError("'trimmed' attribute is not set.")
        data = pd.DataFrame()
        data['cast_id'] = [self.cast_id]
        best_full_cast = self.filtered if self.filtered is not None else self.proc
        # TODO: the below uses all very odf/go-ship/seabird-specific column labels
        #   which should ultimately be replaced with standardized names.
        data['start_time'] = [float(self.trimmed["scan_datetime"].head(1))]
        p_max_ind = best_full_cast[self.p_col].argmax()
        data['bottom_time'] = [float(best_full_cast["scan_datetime"][p_max_ind])]
        data['end_time'] = [float(best_full_cast["scan_datetime"].tail(1))]
        data['start_pressure'] = [float(np.around(self.trimmed[self.p_col].head(1), 4))]
        data['max_pressure'] = [float(np.around(self.trimmed[self.p_col].max(), 4))]
        if 'ALT' in self.proc.columns:
            data['altimeter_bottom'] = [float(np.around(best_full_cast["ALT"][p_max_ind], 4))]
        if all(col in self.proc.columns for col in ['GPSLAT', 'GPSLON']):
            data['latitude'] = [float(np.around(best_full_cast["GPSLAT"][p_max_ind], 4))]
            data['longitude'] = [float(np.around(best_full_cast["GPSLON"][p_max_ind], 4))]
        return data

    def get_pressure_offsets(self, data, threshold, sample_freq):
        """
        Get starting and ending on-deck pressure averages for the cast. A
        conductivity threshold is used to determine when the instrument
        is in or out of the water. A default conductivity threshold is
        set in the user configuration file.

        Previously was 'remove_on_deck' in ctdcal.process_ctd.

        Parameters
        ----------
        data : DataFrame
            Cast data.
        threshold : float
            Maximum conductivity value to determine when cast is not in the water,
            in source conductivity units.
        sample_freq : int
            Instrument sample frequency in Hz

        Returns
        -------
        DataFrame
        """
        # TODO: review algorithm for soundness and reliability, see if we can
        #     use pressure or something other than salinity, remove reliance
        #     on WOCE colnames

        # Half minute
        time_delay = sample_freq * 30  # time to let CTD pressure reading settle/sit on deck
        # split dataframe into upcast/downcast
        downcast = data.iloc[: (data[self.p_col].argmax() + 1)]
        upcast = data.iloc[(data[self.p_col].argmax() + 1):]
        # Search each half of df for minimum conductivity
        # threshold to identify when rosette is out of water
        start_df = downcast.loc[
            (downcast['CTDCOND1'] < threshold)
            & (downcast['CTDCOND2'] < threshold),
            self.p_col,
        ]
        end_df = upcast.loc[
            (upcast['CTDCOND1'] < threshold)
            & (upcast['CTDCOND2'] < threshold),
            self.p_col,
        ]
        # Evaluate starting and ending pressures
        start_samples = len(start_df)
        if start_samples > time_delay:
            start_p = np.average(start_df.iloc[(sample_freq * 2): (start_samples - time_delay)])
        else:
            start_seconds = start_samples / sample_freq
            log.warning(
                    f"{self.cast_id}: Only {start_seconds:0.1f} seconds of start pressure averaged."
            )
            start_p = np.average(start_df.iloc[(sample_freq * 2):start_samples])
        end_samples = len(end_df)
        if end_samples > time_delay:
            end_p = np.average(end_df.iloc[time_delay:])
        else:
            end_seconds = end_samples / sample_freq
            log.warning(
                    f"{self.cast_id}: Only {end_seconds:0.1f} seconds of end pressure averaged."
            )
            end_p = np.average(end_df)  # just average whatever there is
        # Warn if failures
        if len(start_df) == 0:
            log.warning("Failed to find starting deck pressure.")
        if len(end_df) == 0:
            log.warning("Failed to find ending deck pressure.")

        df = pd.DataFrame()
        df['cast_id'] = [self.cast_id]
        df['pressure_start'] = [start_p]
        df['pressure_end'] = [end_p]
        return df

    def trim_ondeck(self, data, threshold):
        """
        Trims end of casts based on a conductivity threshold. Legacy functionality.

        Parameters
        ----------
        data : DataFrame
            Cast data.
        threshold : int
            Maximum conductivity value to determine when cast is not in the water,
            in source conductivity units.
        """
        # split dataframe into upcast/downcast
        # split dataframe into upcast/downcast
        downcast = data.iloc[: (data[self.p_col].argmax() + 1)]
        upcast = data.iloc[(data[self.p_col].argmax() + 1):]
        # Search each half of df for minimum conductivity
        # threshold to identify when rosette is out of water
        start_df = downcast.loc[
            (downcast['CTDCOND1'] < threshold)
            & (downcast['CTDCOND2'] < threshold),
            self.p_col,
        ]
        end_df = upcast.loc[
            (upcast['CTDCOND1'] < threshold)
            & (upcast['CTDCOND2'] < threshold),
            self.p_col,
        ]
        self.ondeck_trimmed = data.iloc[start_df.index.max() : end_df.index.min()].copy()
