import csv
import pandas as pd
import numpy as np
import sys
sys.path.append('ctdcal/')
import settings

def reft_loader(ssscc, reft_dir):

    reft_path = reft_dir + ssscc + '.cap'
    with open(reft_path, 'r', newline='') as f:
        reftF = csv.reader(f, delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
        
        # read in data and convert to dataframe
        reftArray = []
        for row in reftF:
            if len(row) != 17: # skip over 'bad' rows (empty lines, comments, etc.)
                continue
            reftArray.append(row)
        reftArray = np.array(reftArray)
        reftDF = pd.DataFrame(reftArray)
        
        # remove text columns, only need numbers and dates
        reftDF.replace(to_replace=['bn','diff','val','t90','='], value=np.nan, inplace=True)        
        reftDF.dropna(axis=1,inplace=True)

        # combine date/time columns (day/month/year/time are read in separately)
        reftDF[1] = reftDF[1] + ' ' + reftDF[2] + ' ' + reftDF[3] + ' ' + reftDF[4]
        reftDF.drop(columns=[2,3,4], inplace=True)
        
        # rename columns and recast datatypes
        col_names = ['index_memory', 'datetime', 'btl_fire_num', 'diff', 'raw_value', 'T90']
        reftDF.columns = col_names
        reftDF = reftDF.astype({'index_memory':np.int32, 'datetime':object, \
                    'btl_fire_num':np.int32, 'diff':np.int32, \
                    'raw_value':np.float64, 'T90':np.float64})

        # assign initial qality flags
        reftDF.loc[:,'REFTMP_FLAG_W'] = 2
        reftDF.loc[abs(reftDF['diff']) >= 3000, "REFTMP_FLAG_W"] = 3

        # add in STNNBR, CASTNO columns
        reftDF['STNNBR'] = ssscc[0:3]
        reftDF['CASTNO'] = ssscc[3:5]

    return reftDF

def process_reft(ssscc, reft_dir):

    try: 
        reftDF = reft_loader(ssscc, reft_dir)
        reftDF.to_csv(reft_dir + ssscc + '_reft.csv', index=False)
    except FileNotFoundError:
        print('refT file for cast ' + ssscc + ' does not exist... skipping')
        return
