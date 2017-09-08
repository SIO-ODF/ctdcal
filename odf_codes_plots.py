import matplotlib
matplotlib.use('ps')

import libODF_merge_codes as merge
import libODF_ctd_plots as ctd_plots

import sys
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt

# std_dir = f'/Users/jgum/work_code/cruises/'
# whp_btl_name = f'320620170820_hy1.csv'
#
# cruise_dir = f'{std_dir}NBP1707/'
# qual_codes_filepath = f'{cruise_dir}quality_codes/'
# log_dir = f'{cruise_dir}ctd_proc_rewrite/data/logs/quality_code/all/'
#
# file_ssscc = f'{cruise_dir}ssscc.csv'
# bottle_file= f'{qual_codes_filepath}{whp_btl_name}'

def main(argv):
    std_dir = f'/Users/jgum/work_code/cruises/'
    whp_btl_name = f'320620170820_hy1.csv'

    cruise_dir = f'{std_dir}NBP1707/'
    qual_codes_filepath = f'{cruise_dir}quality_codes/'
    log_dir = f'{cruise_dir}ctd_proc_rewrite/data/logs/quality_code/all/'

    file_ssscc = f'{cruise_dir}ssscc.csv'
    bottle_file= f'{qual_codes_filepath}{whp_btl_name}'

    ### compile bottle trips into one file, then merge with odf_db bottle file
    df_bottle = merge.prelim_ctd_bottle_df(file_ssscc, bottle_file)
    ### compute residuals for plotting
    merge.ctd_residuals_df(df_bottle)
    df_bottle = merge.merge_ctd_codes(log_dir, df_bottle)
    df_coded = merge.get_qc_all_from_df(df_bottle)

    #turn into method later
    df_coded = merge.discount_test_stations(df_coded)
    ctd_plots.all_plots(df_coded)

if __name__ == '__main__':
    main(sys.argv[1:])
