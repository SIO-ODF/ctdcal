import pandas as pd
import ctdcal.merge_codes as merge_codes

std_dir = f'/Users/jgum/work_code/ctdcal/'
whp_btl_name = f'320620170820_hy1.csv'

cruise_dir = f'{std_dir}NBP1707/'
qual_codes_filepath = f'{cruise_dir}quality_codes/'
log_dir = f'{cruise_dir}ctd_proc_rewrite/data/logs/quality_code/all/'

file_ssscc = f'{cruise_dir}ssscc.csv'
bottle_file= f'{qual_codes_filepath}{whp_btl_name}'

def write_df_pkl(df, filename):
    df.to_pickle(f'{cruise_dir}ctd_proc_rewrite/data/bottle/{filename}.pkl')
    return None

def main(argv):
    df_bottle_trip = merge_codes.merge_bottle_trip_dfs(file_ssscc)

    write_df_pkl(df_bottle_trip, 'bottle_raw')
    write_df_pkl(df_bottle_trip, 'bottle_adjusted')

if __name__ == '__main__':
    main(sys.argv[1:])
