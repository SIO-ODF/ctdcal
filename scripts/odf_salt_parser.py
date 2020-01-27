import csv
import pandas as pd
import numpy as np
import ctdcal.fit_ctd as fit_ctd
import os

def salt_loader(saltpath):
    f = open(saltpath, newline='')
    saltF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
    
    saltArray = []
    for row in saltF:
        saltArray.append(row)
    del saltArray[0]
         
    header = ['STNNBR','CASTNO','SAMPNO','BathTEMP','CRavg','autosalSAMPNO',\
              'Unknown','StartTime','EndTime','Attempts','Reading1','Reading2',\
              'Reading3', 'Reading4', 'Reading5','Reading6','Reading7','Reading8',\
              'Reading9', 'Reading10','Reading11','Reading12']
    f.close()
    # make all rows of Salt files the same length as header   
    for row in saltArray:
        if len(row) < len(header):
            row.extend([np.NaN]*(len(header)-len(row)))
            
    saltArray = np.array(saltArray) # change to np array
    
    saltDF = pd.DataFrame(saltArray,columns=header) # change to DataFrame
    saltDF = saltDF.apply(pd.to_numeric, errors='ignore')
    saltDF.replace(to_replace='nan', value=np.NaN,inplace=True)
    saltDF.dropna(axis=1,how='all',inplace=True)
    saltDF = saltDF[saltDF['autosalSAMPNO'] != 'worm']
    saltDF['STNNBR'] = saltDF['STNNBR'].astype(int) # force to int (if loaded as str)
    saltDF['SALNTY'] = fit_ctd.SP_salinometer((saltDF['CRavg']/2.0),saltDF['BathTEMP'])
    return saltDF

def salt_df_parser(saltDF, outdir, stn_col='STNNBR', cast_col='CASTNO'):
    stations = saltDF[stn_col].unique()
    for station in stations:
        saltStation = saltDF[saltDF[stn_col] == station]
        casts = saltStation[cast_col].unique()
        for cast in casts:
            stn_cast_salts = saltStation[saltStation[cast_col] == cast]
            # format to SSSCC_salts.csv
            outfile = outdir + '{0:03}'.format(station) + '{0:02}'.format(cast) + '_salts.csv'
            if not os.path.exists(outfile):
                stn_cast_salts.to_csv(outfile,index=False)
            else:
                print(outfile + ' already exists...skipping')