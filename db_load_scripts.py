'''Scripts to load csv data into db.
Expand to numpy arrays/pandas dataframes later at standardization step.



'''
import db_lib
import sys
from datetime import datetime
import csv
from pytz import timezone


def load_BottleFire(f):
    '''Load bottle fire location info. Station/cast/position/time/pressure only.
    '''
    

def load_SBE35(f, station, cast):
    '''Compares time against BottleFire table to find station/cast info, then upload to db.
    Must be run after corresponding bottle file has been uploaded to db, or will fail.
    '''
    session = db_lib.start_talk()

    #Transpose vectors. When written in numpy arrays or pandas this will be easier
    a1 = []
    a2 = []
    a3 = []
    a4 = []
    a5 = []
    a6 = []
    a_station = []
    a_cast = []
    csvr = csv.reader(f)
    next(csvr, None)
    for row in csvr:
        for i, x in enumerate(row):
            if i == 0:
                a1.append(x)
                a_station.append(station)
                a_cast.append(cast)
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
    #for i, j, k ,l in zip(a_station, a_cast, a3, a6):
    #    print(i, j, k, l)
    '''
    #attempt at auto determining with time. Timestamps between bottle file and sbe35
    #are too different to be able to determine, so the user must manually input
    #station cast each time. Revisit this code in the future possibly.

    a_pressure = []

    for i, x in enumerate(a2):
        #Example for format: '20 Feb 2016 10:46:28'
        a2[i] = datetime.strptime(x,'%d %b %Y %H:%M:%S').replace(tzinfo=timezone('UTC')).timestamp()

    for x, y in zip(a2, a3):
        a_pressure.append(db_lib.getSBE35pressure(session, station, cast, y))
    for x, y in zip(a2, a3):
        print(x)
    '''
    db_lib.setSBE35Bulk(session, a_station, a_cast, a3, a6)
    db_lib.end_talk()
    return True

def cut_SBE35(f, fname):
    '''Cuts a SBE35 file into multiple files with properly named files.
    It assumes the first file name is named with the correct station/cast, and then
    queries for the next station/cast that is a CTD cast, and names it accordingly.
    '''
    #TODOLATER
    return None
