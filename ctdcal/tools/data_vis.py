"""
A quick visualizer of just the data presented in bottle file of the processing folder.
"""

#   1. Initialize
import pandas as pd
import numpy as np
import glob

#   2. Setup data matrices
# load bottle file (by MK)
fname = glob.glob("../data/pressure/*hy1.csv")[0]    #   fname is retrning the filename of the (first) bottle file (not use ../.. if running from tools\)
btl_data = pd.read_csv(fname, skiprows=[0, 2], skipfooter=1, engine="python", na_values=-999)   #   read the csv as btl_data. Rows 0 and 2 are headers/units. Footer is blank.
btl_data["SSSCC"] = btl_data["STNNBR"].apply(lambda x: f"{x:03d}") + btl_data[
    "CASTNO"
].apply(lambda x: f"{x:02d}")   #   Create an SSSCC column, assembled from the STNNBR and CASTNO. The lambda functions are to append zeros where missing (1 becomes 001 for SSS, 1 becomes 01 for CC)
#   If you were to just add the strings together, you'll get a new line character for each entry.


btl_data["Residual"] = btl_data["SALNTY"] - btl_data["CTDSAL"]
btl_data[["CTDPRS", "Residual"]] = btl_data[["CTDPRS", "Residual"]].round(4)    #   Round to 4 decimal points
#   Should probably use numpy to do any data touchups/stats

print(btl_data.head())
#   3. Setup Bokeh

#   Cumulative cross section
#   Geographic
#   Cast property property
#   Table of values


#   4. Run Bokeh