#!/usr/bin/env python
import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def report_pressure_details(stacast, log_file, start, end):
    """
    Write start/end deck pressure to ondeck_pressure.csv log file.

    Parameters
    ----------
    stacast : str
        station cast data for file
    c_file : str
        file name location for cast details
    start : str
        cast start time from top of cast after 10m soak
    end : str
        cast end time when instrument leaves water

    Returns
    -------
    None
    """
    df = pd.DataFrame(
        {"SSSCC": stacast, "ondeck_start_p": start, "ondeck_end_p": end}, index=[0]
    )
    add_header = not Path(log_file).exists()  # add header iff file doesn't exist
    with open(log_file, "a") as f:
        df.to_csv(f, mode="a", header=add_header, index=False)

    return True


def report_cast_details(
    stacast, c_file, start, end, bottom, start_p, max_p, b_alt, b_lat, b_lon
):
    """
    Write cast details to cast_details.csv log file.

    Parameters
    ----------
    stacast : str
        station cast data for file
    c_file : str
        file name location for cast details
    start : str
        cast start time from top of cast after 10m soak
    end : str
        cast end time when instrument leaves water
    bottom : str
        bottom of cast time when instrument reaches max depth
    start_p : str
        starting pressure at the time the cast begins
    max_p : str
        maximum pressure for entire cast
    b_alt : str
        altimeter value at bottom of cast
    b_lat : str
        latitude at bottom of cast
    b_lon : str
        longitude at bottom of cast

    Returns
    -------
    None
    """
    df = pd.DataFrame(
        {
            "SSSCC": stacast,
            "start_time": start,
            "bottom_time": bottom,
            "end_time": end,
            "start_pressure": start_p,
            "max_pressure": max_p,
            "altimeter_bottom": b_alt,
            "latitude": b_lat,
            "longitude": b_lon,
        },
        index=[0],
    )
    add_header = not Path(c_file).exists()  # add header iff file doesn't exist
    with open(c_file, "a") as f:
        df.to_csv(f, mode="a", header=add_header, index=False)

    return
