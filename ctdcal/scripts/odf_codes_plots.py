import matplotlib

from . import ctd_plots as ctd_plots
from . import merge_codes as merge

matplotlib.use("ps")

expocode = "320620180309"

qual_codes_filepath = "data/quality_codes/"
log_dir = "data/logs/quality_code/all/"

file_ssscc = "data/ssscc.csv"
bottle_file = f"{qual_codes_filepath}{expocode}_hy1.csv"

### compile bottle trips into one file, then merge with odf_db bottle file
df_bottle = merge.prelim_ctd_bottle_df(file_ssscc, bottle_file)
### compute residuals for plotting
merge.ctd_residuals_df(df_bottle)
df_bottle = merge.merge_ctd_codes(log_dir, df_bottle)
df_coded = merge.get_qc_all_from_df(df_bottle)
print(merge.residual_stddev(df_coded))

##### PATCHED IN FOR NBP1802 #####
df_coded = df_coded.rename(index=str, columns={"CTDCOND2_FLAG_W": "CTDSAL_FLAG_W"})
df_coded.to_pickle("data/values_and_codes.pkl")
# print(df_coded.axes)
##### PATCHED IN FOR NBP1802 #####

ctd_plots.make_cruise_report_plots(df_coded)
