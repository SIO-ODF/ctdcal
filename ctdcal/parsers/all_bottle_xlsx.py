"""
A quick & dirty parser to make bottle sample csv files from an excel spreadsheet.
Developed for use during ICES intercomparison project.
"""
from pathlib import Path

import pandas as pd

from ctdcal.common import validate_dir

def parse_discrete(infile, cnvdir, name, cast_list, colname, export_colname, cast_id_col="Cast", btlnum_col="Bottle Number"):
    outdir = validate_dir(Path(cnvdir), create=True)
    data = pd.read_excel(infile, usecols=[cast_id_col, btlnum_col, colname], dtype={cast_id_col: str})
    data.rename(columns={btlnum_col: "SAMPNO", colname: export_colname}, inplace=True)
    for cast_id in cast_list:
        outfile = Path(outdir, '%s_%s.csv' % (cast_id, name))
        outdata = data.loc[(data[cast_id_col] == cast_id) & (data[export_colname].notnull())]
        if len(outdata) > 0:
            outdata.to_csv(outfile, columns=["SAMPNO", export_colname], index=False)
