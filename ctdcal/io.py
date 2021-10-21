from pathlib import Path
from typing import Union

import pandas as pd


def load_exchange_btl(btl_file: Union[str, Path]) -> pd.DataFrame:
    """
    Load WHP-exchange bottle file (_hy1.csv) into DataFrame.

    Parameters
    ----------
    btl_file : str or Path
        Name of file to be loaded

    Returns
    -------
    df : DataFrame
        Loaded bottle file
    """
    with open(btl_file) as f:
        file = f.readlines()
        for idx, line in enumerate(file):
            if line.startswith("EXPOCODE"):
                units = idx + 1  # units row immediately follows column names

    return pd.read_csv(
        btl_file,
        skiprows=[0, units],
        skipfooter=1,
        engine="python",
        comment="#",
        skipinitialspace=True,
    )


def write_pressure_details(
    ssscc: str, log_file: Union[str, Path], p_start: float, p_end: float
) -> None:
    """
    Write start/end deck pressure to ondeck_pressure.csv log file.

    Parameters
    ----------
    ssscc : str
        Station/cast in SSSCC format
    log_file : str or Path
        File destination for pressure details
    p_start : float
        Average starting on-deck pressure (pre-deployment)
    p_end : float
        Average ending on-deck pressure (post-deployment)

    Returns
    -------
    None
    """
    df = pd.DataFrame(
        {"SSSCC": ssscc, "ondeck_start_p": p_start, "ondeck_end_p": p_end}, index=[0]
    )
    add_header = not Path(log_file).exists()  # add header iff file doesn't exist
    with open(log_file, "a") as f:
        df.to_csv(f, mode="a", header=add_header, index=False)


def write_cast_details(
    ssscc: str,
    log_file: Union[str, Path],
    time_start: float,
    time_end: float,
    time_bottom: float,
    p_start: float,
    p_max: float,
    b_alt: float,
    b_lat: float,
    b_lon: float,
) -> None:
    """
    Write cast details to cast_details.csv log file.

    Parameters
    ----------
    ssscc : str
        Station/cast in SSSCC format
    log_file : str or Path
        File destination for cast details
    time_start : float
        Time at start of cast (from minimum pressure after 10m soak)
    time_end : float
        Time at end of cast (when instrument leaves water)
    time_bottom : float
        Time at bottom of cast (max depth)
    p_start : float
        Pressure at the time the cast begins
    p_max : float
        Pressure at bottom of cast
    b_alt : float
        Altimeter value at bottom of cast
    b_lat : float
        Latitude at bottom of cast
    b_lon : float
        Longitude at bottom of cast

    Returns
    -------
    None
    """
    df = pd.DataFrame(
        {
            "SSSCC": ssscc,
            "start_time": time_start,
            "bottom_time": time_bottom,
            "end_time": time_end,
            "start_pressure": p_start,
            "max_pressure": p_max,
            "altimeter_bottom": b_alt,
            "latitude": b_lat,
            "longitude": b_lon,
        },
        index=[0],
    )
    add_header = not Path(log_file).exists()  # add header iff file doesn't exist
    with open(log_file, "a") as f:
        df.to_csv(f, mode="a", header=add_header, index=False)
