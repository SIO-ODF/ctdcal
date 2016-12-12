'''Scripts to load csv data into db.
Expand to numpy arrays/pandas dataframes later at standardization step.



'''
import db_lib
import sys
from datetime import datetime
import csv

def load_BottleFire(f):
    '''Load bottle fire location info. Station/cast/position/time/pressure only.
    '''

def load_SBE35(f):
    '''Compares time against BottleFire table to find station/cast info, then upload to db.
    Must be run after corresponding bottle file has been uploaded to db, or will fail.
    '''
    #Transpose vectors. When written in numpy arrays or pandas this will be easier
    a1 = []
    a2 = []
    a3 = []
    a4 = []
    a5 = []
    a6 = []
    csvr = csv.reader(f)
    next(csvr, None)
    for row in csvr:
        for i, x in enumerate(row):
            if i == 0:
                a1.append(x)
            if i == 1:
                a2.append(x)
            if i == 2:
                a3.append(x)
            if i == 3:
                a4.append(x)
            if i == 4:
                a5.append(x)
            if i == 5:
                a6.append(x)

    #a7 = []
    for i, x in enumerate(a2):
        #20 Feb 2016 10:46:28
        a2[i] = datetime.strptime(x,'%d %b %Y %H:%M:%S')
        #a7.append(datetime.strptime(x,'%d %b %Y %H:%M:%S'))
    for x in a2:
        print(x)

    
