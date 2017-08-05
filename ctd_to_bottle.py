import pandas as pd
import numpy as np
import io
import sys

def main(argv):
    '''Creates a bottle file with CTD downcast information.

    Needs to be cleaned up from script form to other form, hopefully
    '''
    #general fileloading area
    output_file = '/Users/jgum/work_code/cruises/P06_2017/ctd_proc_rewrite/data/scratch_folder/ctd_to_bottle.csv'
    file_ssscc = '/Users/jgum/work_code/cruises/P06_2017/ssscc.csv'
    df_all = pd.DataFrame()

    all_casts = []
    with open(file_ssscc, 'r') as filename:
        all_casts = [line.strip() for line in filename]

    for ssscc in all_casts:
        station = int(ssscc[0:3])
        cast = int(ssscc[3:5])
        #bottle handling section
        dir_bottle = '/Users/jgum/work_code/cruises/P06_2017/ctd_proc_rewrite/data/bottle/'
        bottle_postfix = '_btl_mean.csv'
        btl_skiprows = [0,1]
        btl_names = ['SAMPNO','CTDPRS']
        btl_usecols = [0,3]

        df_bottle = pd.read_csv(dir_bottle + ssscc + bottle_postfix, skiprows=btl_skiprows, names=btl_names, usecols = btl_usecols)
        df_bottle['bins'] = pd.cut(df_bottle['CTDPRS'], range(0,int(np.ceil(df_bottle['CTDPRS'].max()))+5,2), right=False, include_lowest=True)

        #ctd handling section
        dir_ctd = '/Users/jgum/work_code/cruises/P06_2017/ctd_proc_rewrite/data/pressure/'
        ctd_postfix = '_ct1.csv'
        ctd_skiprows = [0,1,2,3,4,5,6,7,8,9,10,11,13]

        df_c = pd.read_csv(dir_ctd + ssscc + ctd_postfix, skiprows=ctd_skiprows, skipfooter=1, engine='python')
        df_c['bins'] = pd.cut(df_c['CTDPRS'],range(0,int(np.floor(df_c['CTDPRS'].max()))+5,2), right=False, include_lowest=True) #heisenbug here at bin padding

        df_merged_btl_ctd_4 = pd.merge(df_c, df_bottle[['bins']], on='bins', how='outer')

        df_merged_btl_ctd_6 = pd.merge(df_merged_btl_ctd_4.iloc[:,1:], df_bottle, on='bins', how='inner')
        df_merged_btl_ctd_6.groupby('SAMPNO').first()

        #reference temperature section
        #there should not be any 9 values, unless the reftemp dies? then we have a PROBLEM
        dir_reft = '/Users/jgum/work_code/cruises/P06_2017/ctd_proc_rewrite/data/reft/'
        reft_postfix = '_reft.csv'
        try:
            df_reft = pd.read_csv(dir_reft + ssscc + reft_postfix, names =['SAMPNO', 'REFTMP'], skiprows = 1, usecols = [2, 5])
            df_reft['REFTMP_FLAG_W'] = 2
            df_incomplete_reft = pd.merge(df_reft, df_merged_btl_ctd_6, on='SAMPNO', how='outer')
            df_incomplete_reft.groupby('SAMPNO').fillna(inplace=True, method = 'pad')
            df_complete_reft = df_incomplete_reft.fillna(value={'REFTMP':-999, 'REFTMP_FLAG_W':5})
        except FileNotFoundError:
            df_merged_btl_ctd_6['REFTMP'] = -999
            df_merged_btl_ctd_6['REFTMP_FLAG_W'] = 5
            df_complete_reft = df_merged_btl_ctd_6

        df_complete_reft['STNNBR'] = station
        df_complete_reft['CASTNO'] = cast
        df_complete_reft = df_complete_reft.groupby('SAMPNO').first()
        df_complete_reft.reset_index(inplace=True)
        df_all = df_all.append(df_complete_reft)

    #create the stupid string
    #d_start = {'0':['CASTNO,CTDOXY,CTDOXY_FLAG_W,CTDPRS,CTDPRS_FLAG_W,CTDSAL,CTDSAL_FLAG_W,CTDTMP,CTDTMP_FLAG_W,CTDFLUOR,CTDFLUOR_FLAG_W,REFTMP,REFTMP_FLAG_W,SAMPNO,STNNBR,CTDXMISS,CTDXMISS_FLAG_W',',UMOL/KG,,DBAR,,PSS-78,,ITS-90,,VOLTS,,ITS-90,,,,VOLTS,']}
    #d_end = {'blah' : ['END_DATA']}
    #df_start = pd.DataFrame(d_start)
    #df_end = pd.DataFrame(d_end)
    stringy = 'CASTNO,CTDOXY,CTDOXY_FLAG_W,CTDPRS,CTDPRS_FLAG_W,CTDSAL,CTDSAL_FLAG_W,CTDTMP,CTDTMP_FLAG_W,CTDFLUOR,CTDFLUOR_FLAG_W,REFTMP,REFTMP_FLAG_W,SAMPNO,STNNBR,CTDXMISS,CTDXMISS_FLAG_W\n,UMOL/KG,,DBAR,,PSS-78,,ITS-90,,VOLTS,,ITS-90,,,,VOLTS,\n'
    #stringy = df_start.to_string(header=False, index=False) + '\n'
    #the next line is incredibly stupid. find a better way to create
    stringy += df_all.iloc[:,0:-1].to_string(header=False, index=False, formatters=[lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)+',',lambda x: str(x)])
    #stringy += df_end.to_string(header=False, index=False)
    stringy += '\nEND_DATA\n'

    with open(output_file, 'w') as f_handle:
        f_handle.write(stringy)

if __name__ == '__main__':
    main(sys.argv[1:])