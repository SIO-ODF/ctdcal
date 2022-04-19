import logging

# import warnings
from datetime import datetime, timezone
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
from py import process
import scipy.signal as sig

from ctdcal import get_ctdcal_config, io, oxy_fitting, process_ctd

import pytest



def writer(tmp_path, stringname, data):
    """A function for writing a temporary testing file on tmp_path,
    called 'stringname', using 'data"""
    with open(tmp_path / stringname, "w+") as f:
        for line in data:
            f.write(line)


def create_df(path, rows=4001, n_deck=10, n_soak=200, p_outliers=0):
    #   Below is what MK did: Try using this as an example to set up a DataFrame for testing
    #   What do you need inside of it? Will it error out?
    #  import numpy as np
    #  # make fake data
    # downcast = np.concatenate(
    #     (
    #         [0] * n_deck,  # on deck
    #         np.linspace(0, 20, n_soak),
    #         [20] * n_soak * 3,  # fake soak
    #         np.linspace(20, 0, n_soak),
    #         np.linspace(0, 500, rows - n_soak * 5),
    #     )
    # )
    # upcast = np.linspace(downcast[-1], 0, n_deck + rows - 1)
    # cond = np.concatenate(([0] * n_deck, np.linspace(30, 35, rows)))
    # temp = np.concatenate(([10] * n_deck, np.linspace(5, 0, rows)))
    # pump = np.concatenate(([0] * n_deck, [1] * rows))

    # # construct and export DataFrame
    # df = pd.DataFrame()
    # df["CTDPRS"] = np.concatenate((downcast, upcast))
    # df["CTDCOND1"] = np.concatenate((cond, cond[:-1][::-1]))
    # df["CTDCOND2"] = np.concatenate((cond, cond[:-1][::-1]))
    # df["CTDSAL"] = np.concatenate((cond, cond[:-1][::-1]))
    # df["CTDTMP1"] = np.concatenate((temp, temp[:-1][::-1]))
    # df["CTDTMP2"] = np.concatenate((temp, temp[:-1][::-1]))
    # n_data = len(df)

    # # add pressure outlier
    # if p_outliers:
    #     bad_rows = np.random.randint(0, n_data, 2)
    #     df.loc[bad_rows, "CTDPRS"] = 8000

    # # extra weirdness to make remove_on_deck not fail
    # # (TODO: remove need for this junk!!!)
    # df["pump_on"] = np.concatenate((pump, pump[:-1][::-1]))
    # df["scan_datetime"] = 1645231099 + np.linspace(0, 1e3, n_data)
    # df["GPSLAT"] = np.array([45] * n_data)
    # df["GPSLON"] = 1645231099 + np.linspace(0, 1e3, n_data)
    # df["ALT"] = 590 - df["CTDPRS"]  # altimeter goes 100m->10m and back
    # df.loc[df["ALT"] > 100, "ALT"] = 100

    # df.to_pickle(path)
    pass

#  Test conditions: Correct csv, one station, one station (no \n), strip() test
@pytest.mark.parametrize(
    "data,filename",
    [
        (["03001\n", "04001\n"], "good1"),
        ("04040\n", "good2"),
        ("12345", "good3"),
        (["  01001\n  ", " 02001 \n", "03001\n", " 04001 \n "], "good4"),
    ],
)
def test_get_ssscc_list(caplog, tmp_path, data, filename):
    """
    * If it can't find the right directory, have it error out gracefully
        * Assert that the errors are what we expect (missing dir)
    * If the directory list is faulty, have it error out gracefully
        * Assert that the file is read in as a list
    """

    #  Run without any sort of file structure and have it return FileNotFoundError
    with pytest.raises(FileNotFoundError):
        ssscc = process_ctd.get_ssscc_list()

    #  Run parametrized test scenarios
    to_write = filename + "_ssscc.csv"
    writer(tmp_path, to_write, data)
    with caplog.at_level(logging.INFO):
        ssscc = process_ctd.get_ssscc_list(tmp_path / to_write)
    assert type(ssscc) == list, "ssscc is a list"
    assert len(ssscc) != 0, "The list is not empty"
    # print(ssscc[0])
    # assert False  #  Should you want to see the stout (print statements)


CONTENT = "content"


def test_make_ssscc_list(tmpdir):
    #  Run without any file structure and have it throw a FileNotFoundError
    with pytest.raises(FileNotFoundError):
        ssscc = process_ctd.make_ssscc_list()

    #  Create /data/raw and write files to it
    a = tmpdir.mkdir("data")  #  Mandatory to make subfolder /raw/
    b = tmpdir.mkdir("data/raw")
    # assert bool(process_ctd.make_ssscc_list(a + '/empty.csv', Path(b))) == 0 #  Confirm empty list returned
    filenames = ["01001.hex", "02001.hex", "03001.hex"]
    for name in filenames:
        x = b.join(name)
        x.write("contents")
    assert x.read() == "contents", "Last file written to correctly."
    
    #   Build the ssscc.csv and test it
    # files = Path(b).glob('*.hex')
    ssscc = process_ctd.make_ssscc_list(a + "/ssscc.csv", Path(b))
    assert type(ssscc) == list  #  Check return
    assert bool(ssscc)
    assert ssscc[0] == filenames[0].replace('.hex', '') #  Should be just the #

def test_trim_soak_period():
    #   Test and raise any expected errors with an empty DF

    #   Define a DataFrame and pass it in

    #   Ensure output DataFrame is smaller than before

    pass

def test_cast_details():
    #   Define a DataFrame
    #   Define SSSCC for the DataFrame

    #   WRITE SEPERATE TESTS FOR:
        #   _trim_soak_period
        #   io.write_cast_details (already written by MK)

    #   Get back df_downcast dataframe:
    #   * Confirm that it's not empty
    #   * Confirm that the data was trimmed (len())
    pass
