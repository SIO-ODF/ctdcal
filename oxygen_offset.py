import pandas as pd
import libODF_merge_codes as merge
import numpy as np
import io

qual_codes_filepath = f'/Users/jgum/work_code/cruises/NBP1707/quality_codes/'
cruise_dir = f'/Users/jgum/work_code/cruises/NBP1707/'
log_dir = f'{cruise_dir}ctd_proc_rewrite/data/logs/quality_code/all/'

file_ssscc = f'{cruise_dir}ssscc.csv'
bottle_file= f'{qual_codes_filepath}320620170820_hy1.csv'

### compile bottle trips into one file, then merge with odf_db bottle file
#df_bottle = merge.prelim_ctd_bottle_df(file_ssscc, bottle_file)
### compute residuals for plotting
#merge.ctd_residuals_df(df_bottle)
df_bottle = merge.load_exchange_bottle_file(bottle_file)
df_bottle = merge.merge_ctd_codes(log_dir, df_bottle)
#df_bottle = merge.get_qc_all_from_df(df_bottle)

df_bottle['CTDOXY'] = df_bottle['CTDOXY'] - 2.2
df_bottle['CTDOXY_FLAG_W'] = df_bottle['CTDOXY_FLAG_W'] + 1
df_bottle = df_bottle.rename(index=str, columns={'CTDCOND1_FLAG_W': 'CTDSAL_FLAG_W'})

df_bottle.to_pickle('data/values_and_codes.pkl')
print(f'Finished pickling bottle file')
all_casts = []
with open(file_ssscc, 'r') as filename:
    all_casts = [line.strip() for line in filename]

for ssscc in all_casts:
    lines = ''
    with open(f'{cruise_dir}ctd_proc_rewrite/data/scratch_folder/pressure/{ssscc}_ct1.csv', 'r') as f:
        for x in range(0,14):
            lines += f.readline()
    df = pd.read_csv(f'{cruise_dir}ctd_proc_rewrite/data/scratch_folder/pressure/{ssscc}_ct1.csv',
    skiprows=[0,1,2,3,4,5,6,7,8,9,10,11,13], skipfooter=1, engine='python')
    df['CTDOXY'] = df['CTDOXY'] - 2.2
    df['CTDOXY_FLAG_W'] = df['CTDOXY_FLAG_W'] + 1
    lines += df.to_string(header=False, index=False, formatters=[lambda x: str(x)+',',
    lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',
    lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',
    lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',
    lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)])
    lines += f'\nEND_DATA\n'
    with open(f'{cruise_dir}ctd_proc_rewrite/data/scratch_folder/pressure/{ssscc}_ct1.csv', 'w') as f:
        f.write(lines)
