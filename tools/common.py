#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd


class Cast(object):
    """
    Base class to store cast details.
    """
    def __init__(self):
        self.raw = pd.read_pickle('../data/converted/%s.pkl' % self.cast_no)
        self.btl = pd.read_pickle('../data/bottle/%s_btl_mean.pkl' % self.cast_no)
        self.data = pd.DataFrame()
        self.coeffs = {}


class Coefficients(object):
    """
    A Coefficients class with two methods to load/save the serialized calibration coefficients for an instrument.
    """
    def __init__(self, coeff_file):
        """
        Initialize the class with the path to coefficients file and an empty dictionary structure for
        the calibration data
        """
        # set the infile name and path
        self.coeff_file = coeff_file
        self.coeffs = {}

    def load_csv(self):
        """
        Obtain the calibration data for this instrument from a csv calibration file.
        """
        with open(self.coeff_file, 'r') as f:
            for line in f:
                line = line.split(',')
                self.coeffs[line[1].strip()] = float(line[2].strip())

    # def save_coeffs(self):
    #     """
    #     Save the calibration data for this instrument to a JSON data file.
    #     """
    #     with open(self.coeff_file, 'w') as f:
    #         jdata = json.dumps(self.coeffs, cls=NumpyEncoder)
    #         f.write(jdata)
