import pandas as pd
import numpy as np

import libODF_odf_io as odf_io
#import libODF_precalibration as precal

std_dir = f'/Users/jgum/work_code/cruises/'
cruise_dir = f'{std_dir}NBP1707/'
data_dir = f'{cruise_dir}ctd_proc_rewrite/data/'
bottle_dir = f'{data_dir}bottle/'
salt_dir = f'{data_dir}salt/'

ssscc_file = '../ssscc.csv'

ssscc = []
with open(ssscc_file, 'r') as filename:
    ssscc = [line.strip() for line in filename]

corrected_salts = pd.DataFrame()
for x in ssscc:
    #Boilerplate
    df = odf_io.load_short_salts_2017(f'{salt_dir}{x}')
    df = odf_io.salts_create_index_time(df)
    p = odf_io.autosal_drift_calibration(df)
    df = odf_io.autosal_drift_fit(df)
    df = odf_io.compute_salinity(df)
    df = odf_io.formatted_salt_file(df)

    #Put together the frames to make one big one
    corrected_salts = pd.concat([corrected_salts, df], axis = 0)

corrected_salts.to_pickle(f'{data_dir}/scratch_folder/all_corrected_salts.pkl')
