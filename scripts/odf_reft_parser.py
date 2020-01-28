import csv
import pandas as pd
import numpy as np
import sys
sys.path.append('ctdcal/')
import settings

def parse(f, f_out):
    #second half of final output
    array = []
    reader = csv.reader(f, delimiter=' ')

    #break row down to a standardized size
    for row in reader:
        row2 = []
        for x in row:
            if x != '':
                row2.append(x)

        #if not a fixed size, assume it's a comment and to be ignored
        if len(row2) != 17:
            continue

        #start to build list to be written out
        row3 = []
        #hardcoded madness because it probably won't change anytime soon.
        row3.append(row2[0])
        #concatenate date/time, convert to ISO8601 later
        row3.append(row2[1] + ' ' + row2[2] + ' ' + row2[3] + ' ' + row2[4])
        row3.append(row2[7])
        row3.append(row2[10])
        row3.append(row2[13])
        row3.append(row2[16])
        array.append(row3)

    nparray = np.array(array)
    df = pd.DataFrame(nparray, columns=['index_memory', 'datetime', 'btl_fire_num', 'diff', 'raw_value', 'T90'])
    df = df.astype({'index_memory':np.int32, 'datetime':object, 'btl_fire_num':np.int32, 'diff':np.int32, 'raw_value':np.float64, 'T90':np.float64})
    # assign initial flags
    df.loc[:,'REFTMP_FLAG_W'] = 2
    df.loc[abs(df['diff']) >= 3000, "REFTMP_FLAG_W"] = 3
    # add in STNNBR, CASTNO columns
    df['STNNBR'] = f_out[0:3]
    df['CASTNO'] = f_out[3:5]
    return df