"""
A quick & dirty parser to make bottle sample csv files from an excel spreadsheet.
Developed for use during ICES intercomparison project.
"""
from pathlib import Path

import pandas as pd

from ctdcal.common import validate_dir


def parse_salinity(infile, datadir, inst, cast_list, cast_id_col="Cast", btlnum_col="Bottle Number", sal_col="Salinity"):
    outdir = validate_dir(Path(datadir, 'converted', inst), create=True)
    data = pd.read_excel(infile, usecols=[cast_id_col, btlnum_col, sal_col])
    data.rename(columns={btlnum_col: "SAMPNO", sal_col: "SALNTY"}, inplace=True)
    for cast_id in cast_list:
        outfile = Path(outdir, '%s_salts.csv' % cast_id)
        data.loc[(data[cast_id_col] == cast_id) & (data['SALNTY'].notnull())].to_csv(outfile, columns=["SAMPNO", "SALNTY"], index=False)
