'''Library to create SBE .btl equivalent files.

Joseph Gum SIO/ODF
Nov 7, 2016
'''

import io
import numpy as np
import sys
import csv
import datetime
import statistics
import converter_scaffolding as cnv
import pandas as pd


BOTTLE_FIRE_COL_NAME = 'btl_fire'

DEBUG = False

def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Retrieve the bottle data from a converted file.
def retrieveBottleDataFromFile(converted_file, debug=False):

    converted_df = cnv.importConvertedFile(converted_file, DEBUG)
    
    return retrieveBottleData(converted_df, debug)


# Retrieve the bottle data from a dataframe created from a converted file.
def retrieveBottleData(converted_df, debug=False):
    if BOTTLE_FIRE_COL_NAME in converted_df.columns:
        return converted_df.loc[converted_df[BOTTLE_FIRE_COL_NAME] == True]
    else:
        debugPrint("Bottle fire column:", BOTTLE_FIRE_COL_NAME, "not found")
    
    return pd.DataFrame() #empty dataframe



def handler(converted_file, config_file=False, debug=False):
    """Wrapper for the whole thing.
    Take in the file, then call other methods on it to process.
    SBE Data Processing manual does not specify how they average their data,
    so this will need some tweaking to make it an exact copy. Otherwise change
    to fit ODF needs.

    Input:
    converted_file: filename
    config_file: filename for config file to set additional options
    debug: debug setting

    Output:
    output: a string of data to be written to file
    """

    with open(converted_file, 'r') as f:
        datareader = csv.reader(f)
        title = datareader.__next__()
        datareader.__next__()
        header = datareader.__next__()

        #determine bottle_fire column
        btl_index = -1
        for counter, title in enumerate(title):
            if title == 'bottle_fire':
                btl_index = counter
        if btl_index == -1:
            errPrint('bottle_fire not found, please check .converted for bottle_fire column')
            sys.exit(1)
        #the mix of data structures is super messy. figure out some way to simplify it
        #also it makes the input to group_scans() very ehhhhh and needs to be cleaned up
        temp_array = []
        for counter, line in enumerate(datareader):
            if line_extract(line, btl_index):
                temp_array.append((counter, line))

        #return temp_array
        #messy input to group_scans()
        output = group_scans(temp_array, header)
        return output
                #aux struct to hold counter values from group_scans?


def line_extract(row, index):
    """Given a numpy.array row, determine if a bottle has been fired recently.
    Return True or False to tell the upper level system whether to save the row.

    Input:
    row: a numpy.array row equivalent to a single scan of data
    index: an integer for the column to look in for bottle fire information

    Output:
    True/False: boolean whether a row should be saved or not
    """

    try:
        if row[index] == 'True':
            return True
        elif row[index] == 'False':
            return False
    except TypeError:
        debugPrint('Not a boolean - check the column passed to line_extract()')


def group_scans(array, header):
    """Groups up scans by scan count to be passed to bottle_avg().
    Assumes all scans are sequential when looking for values.
    group_scans assumes you want to use only the data with fire bit enabled.
    If you want to use a custom time window use group_scans_custom instead.

    Input:
    array: python list of tuples as (scan, [scan reading])
    header: custom utype superset of numpy to determine how to average readings

    Output:
    output: a string of each bottle data averaged individually to be written to file

    """
    output = ''
    scan = 0
    group_start = 0
    for counter, line in enumerate(array):
        if scan is 0:
            scan = line[0]
            group_start = counter
        elif scan is not 0:
            if line[0] - scan > 1:
                #debugPrint(group_start, counter-1)
                output += bottle_avg(array[group_start:counter-1], header) +'\n'
                group_start = counter
            if line[0] - scan == 1:
                scan = line[0]
            #check this, needs to be fixed later on
            else:
                scan = line[0]
                #print(str(scan) + ': Check input file for errors, scans not sequential')
    return output

def group_scans_custom(array, header, scans):
    """Custom time window version of group_scans.
    Use when a custom time window is wanted, and not the default 1.5-2 second window.
    Assumes the middle of the firing window is the actual bottle close, then computes
    needed scans and reaches forwards and backwards to provide that to bottle_avg.

    Input:
    array: python list of tuples as (scan, [scan reading])
    header: list of custom utype superset of numpy to determine how to average readings
    scans: number of scans total to average, assumes you're using 24Hz

    Output:
    output: a string of data to be written to file

    See group_scans for more information on how group_scans works.

    """
    ##### TO DO #####
    return None


def bottle_avg(array, header):
    """
    Because we're averaging non-numerical data, we'll need handling code for them.
    Ex: Position, datetime, newfix.
    The newfix boolean is treated as an OR, if any boolean is True the avg will be true.

    Input:
    array: array of X lines passed in by group_scans
    header: list of custom utype superset of numpy to determine how to average readings

    Output:
    output: a single string of averaged values
    """
    #create temp list holding all values of a column
    z_list = []
    for counter, row in enumerate(array):
        for inner_counter, value in enumerate(row[1]):
            if counter == 0:
                #print(inner_counter,value)
                z_list.append([value])
            #should be if counter > 0 for clarity
            else:
                #print(z_list)
                #print(counter)
                #print(inner_counter,value)
                z_list[inner_counter].append(value)

    #print(z_list)

    #average values
    temp_out_list = []
    for values, utype in zip(z_list, header):
        temp_out_list.append(average(values, utype))

    #construct output
    output = ''
    for x in temp_out_list:
        output += str(x) + ','
    return output.rstrip(',')

def average(column, utype):
    """
    Overloaded method to handle averaging different types.
    """

    if utype == 'float64':
        temp_list = []
        for x in column:
            temp_list.append(float(x))
        return statistics.mean(temp_list)
    elif utype == 'bool_':
        for x in column:
            if x == 'True':
                return 'True'
        return 'False'
    ###FINISH DATETIME LATER, RETURNS FIRST VALUE FOR NOW
    elif utype == 'datetime':
        return column[0]
    else:
        return None
