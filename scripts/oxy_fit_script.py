#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:36:32 2017

@author: k3jackson
"""
import sys
sys.path.append('../ctdcal/')
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import process_ctd
import pandas as pd
import oxy_fitting
#Make this key argument
method = 3




# These can be input automatically from configuration file

raw_dir = '../data/raw/'
ssscc_file = '../data/ssscc.csv'
time_dir = '../data/time/'
btl_dir = '../data/bottle/'
log_file = '../data/logs/oxy_fit_coefs.csv'



with open(ssscc_file, 'r') as filename:
    ssscc = [line.strip() for line in filename]

dataframe_concat = pd.DataFrame()
coef_concat = pd.DataFrame()

for cast in range(len(ssscc)):

    stn_nbr = int(ssscc[cast][0:3])
    cst_nbr = int(ssscc[cast][-2:])
    stn_cst = str(ssscc[cast])

#    stn_cst = str(ssscc[-3:]+'01')

    hexfile = raw_dir+stn_cst+'.hex'
    xmlfile = raw_dir+stn_cst+'.XMLCON'
    time_file = time_dir+stn_cst+'_time.pkl'
    btl_file = btl_dir+stn_cst+'_btl_mean.pkl'

    time_file = time_dir+stn_cst+'_time.pkl'
    btl_file = btl_dir+stn_cst+'_btl_mean.pkl'

    time_data = process_ctd.dataToNDarray(time_file,float,True,',',1)
    time_data = pd.DataFrame.from_records(time_data)

    btl_data = process_ctd.dataToNDarray(btl_file,float,True,',',0)
    btl_data = pd.DataFrame.from_records(btl_data)


    btl_data_write, btl_data_fit, coef = oxy_fitting.oxy_fit(time_data,btl_data,stn_cst,hexfile,xmlfile,method=method)
    print('COMPLETED Fitting for SSSCC:',stn_cst)
    
    print('Saving Coefficients...')
    print('Saved Coef: ',coef)
    oxy_coef_df = oxy_fitting.write_oxy_coef(coef,stn_cst)
    
    coef_concat = pd.concat([coef_concat,oxy_coef_df])
    dataframe_concat = pd.concat([dataframe_concat,btl_data_fit])

# Save coef to csv

coef_concat.to_csv(log_file)

#df=dataframe_concat
#df['BTL_O'] = df['OXYGEN']-df['CTDOXY']
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#cm = ax.scatter(df['BTL_O'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
#ax.set_xlim(-10,10)
#ax.set_title('OXYGEN-CTDOXY vs CTDPRS')
#ax.set_xlabel('CTDOXY Residual (umol/kg)')
#ax.set_ylabel('Pressure (dbar)')
#cbar = fig.colorbar(cm)
#cbar.set_label('Station Number')
#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#cm = ax.scatter(df['STNNBR'],df['BTL_O'], marker='+', c=df['CTDPRS'], cmap='rainbow')
#ax.set_ylim(-10,10)
#ax.set_title('OXYGEN-CTDOXY vs STNNBR')
#ax.set_xlabel('Station Number')
#ax.set_ylabel('CTDOXY Residual (umol/kg)')
#cbar = fig.colorbar(cm)
#cbar.set_label('Pressure (dbar)')