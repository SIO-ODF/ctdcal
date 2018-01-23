#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 15:10:39 2018

@author: k3jackson
"""
import pandas as pd


def pressure_offset(ondeck_file):
    """
    Needs ondeck pressure csv file
    """
    ondeck_press = pd.read_csv('./ondeck_pressure.csv',header=None,names=['STN','Start','End'])
    ondeck_press['p_start'] = ondeck_press['Start'].str.split(':').str.get(1)
    ondeck_press['p_start'] = ondeck_press['p_start'].apply(pd.to_numeric) 
    ondeck_press['p_end'] = ondeck_press['End'].str.split(':').str.get(1)
    ondeck_press['p_end'] = ondeck_press['p_end'].apply(pd.to_numeric) 
    
    
    return