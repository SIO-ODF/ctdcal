'''To be folded into libODF_process_ctd.py later.'''

import numpy as np
import scipy.signal as sig
import scipy.stats as st
import pandas as pd
import math
import csv #?

def odf_polyfit_handler(master_array, bottle_file, type, degree='auto'):
    '''in case it's not a ndarray, other things'''
    try:

    if type == 'T':
        print('t')
        bottle_array = np.loadtxt(bottle_file, delimiter=',', skiprows=1, usecols = (2,5))

    if type == 'C':
        print('c')


def odf_polyfit(x, y, deg=1):
    '''Polyfit only.'''
    np.polyfit(x,y,deg)
