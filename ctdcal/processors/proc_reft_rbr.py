"""
Processes RBR duet PT sensor data and formats as csv suitable for use
with ctdcal fitting routines.
"""
import logging
from pathlib import Path
from datetime import timedelta as td

import numpy as np
import pandas as pd
from scipy.stats import linregress

from ctdcal.common import validate_dir
from ctdcal.flagging.flag_common import quality_by_threshold
from ctdcal.parsers.parse_rbr_duet_pt import Parser, parse_cast

log = logging.getLogger(__name__)


def proc_reft(cast_list, rawdir, cnvdir, btldir, outdir, export_parsed=False, parsed_dir=None, sync_times=False):
    """
    Processes RBR duet PT sensor data into formatted reference temperature data
    suitable for use with ctdcal fitting routines.

    Parameters
    ----------
    cast_list : list[str]
        list of cast names to process.
    rawdir : str or Path-like
        raw RBR sensor data directory.
    cnvdir : str or Path-like
        processed CTD data directory.
    btldir : str or Path-like
        processed CTD bottle file directory.
    outdir : str or Path-like
        output directory.
    export_parsed : bool
        (optional) if True, save parsed raw RBR sensor data as csv.
    parsed_dir : str or Path-like
        (optional) parsed raw RBR sensor data directory. Required if export_parsed is True.
    """
    if export_parsed is True and parsed_dir is None:
        log.warning("No output directory is set for parsed RBR reft files. Exporting is disabled.")
        export_parsed = False
    elif export_parsed is True:
        parsed_dir = validate_dir(parsed_dir, create=True)
    outdir = validate_dir(outdir, create=True)

    # set up the parser object for parsing the raw files
    parser = Parser(cast_list, rawdir, cnvdir, btldir)

    # loop through each cast with bottle data
    for cast_id in cast_list:
        log.info('Parsing RBR duet PT for cast %s...' % cast_id)
        # parse the rbr data for the cast
        cast = parse_cast(cast_id, parser)
        if cast is None:
            continue
        if export_parsed is True:
            outfile = Path(parsed_dir, '%s.csv' % cast_id)
            cast.data.to_csv(outfile, index=False)

        # apply time synchronization
        if sync_times is True:
            cast.synchronize_cast_data(cnvdir)

        # apply clock drift correction
        cast.data['tstamp'] = pd.to_datetime(
                correct_clock_drift(
                    pd.to_numeric(cast.data['tstamp']).astype('float') / 1E6, parser.meta[cast.raw_file.stem]
                ), unit='ms',
        )

        # get the bottle means for 'interval' seconds after bottle fire times
        parsed_cast_data = cast.data.set_index('tstamp')
        btl_data = get_bottle_data(parsed_cast_data, cast.bottle_times, interval=8, cast_id=cast_id)
        btl_data['cast_id'] = cast_id

        # set some initial quality flags based on std_dev column
        btl_data['quality_flag'] = np.where(
                quality_by_threshold(btl_data['std_dev'], threshold=1),
                3,  # meets or exceeds threshold
                2   # below threshold
        )

        # select and rename columns for output
        output_cols = ['cast_id', 'bottle', 'temp', 'quality_flag']
        processed = btl_data[output_cols].rename(
                columns={'bottle': 'btl_fire_num', 'temp': 'REFTMP', 'quality_flag': 'REFTMP_FLAG_W'}
        )

        # export as csv
        outfile = Path(outdir, '%s_reft.csv' % cast_id)
        processed.to_csv(outfile, index=False)


def get_bottle_data(cast, bottle_times, interval=2, cast_id=None):
    """
    Selects data from cast at bottle fire times and averages over the period
    set by interval. Adds bottle numbers and standard deviation for temperature
    data.
    TODO: needs some data validation e.g. testing bottle times against pressures
        to confirm alignment and correct instrument time (not set to local or
        something), and which pressure to use (pressure_cor vs sea_pressure)

    Parameters
    ----------
    cast : DataFrame
        parsed raw sensor data.
    bottle_times : array
        bottle fire times.
    interval : numeric
        time interval in seconds

    Returns
    -------
    DataFrame
    """
    btl_data = []
    for i, btl_time in enumerate(bottle_times):
        btl_time_start = pd.to_datetime(btl_time, unit='s')
        # set start and stop time
        btl_time_stop = btl_time_start + td(seconds=interval)
        # average the data over the interval period
        btl_interval = cast[btl_time_start: btl_time_stop].mean()[
            ['temp_raw', 'pressure_cor', 'sea_pressure']].to_list()
        # also get a std dev for the temperature col (and pressure)
        btl_interval.append(cast[btl_time_start: btl_time_stop].std()['temp_raw'] * 100)
        btl_interval.append(cast[btl_time_start: btl_time_stop].std()['sea_pressure'] * 100)
        # add the bottle number
        btl_interval.append(i + 1)
        btl_data.append(btl_interval)
    #
    # # # export this untrimmed for troubleshooting alignment (2025-04-21)
    # outdir = '/Users/als026/data/i09n_2025/cruise_data/cnv/rbr'
    # peekfile = Path(outdir, '%s_btl_times.csv' % cast_id)
    # pd.DataFrame(pd.to_datetime(bottle_times, unit='s')).to_csv(peekfile)
    #
    return pd.DataFrame(btl_data, columns=['temp', 'pressure', 'sea_pressure', 'std_dev', 'p_std_dev', 'bottle']).astype(
            {'bottle': int})

def correct_clock_drift(timestamps, meta):
    x = np.array(meta['t_range'])
    y = np.array([meta['t_range'][0], meta['t_range'][1] + meta['t_drift']])
    coeffs = linregress(x, y)
    return coeffs.slope * timestamps + coeffs.intercept
