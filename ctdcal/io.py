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
        btl_file, skiprows=[0, units], skipfooter=1, engine="python", comment="#"
    )
