#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:36:32 2017

@author: k3jackson
"""
import sys
sys.path.append('/Users/k3jackson/p06e/ctd_proc')
sys.path.append('/Users/k3jackson/Kennycode/')
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import libODF_process_ctd as process_ctd
import pandas as pd
import oxy_fitting_ver2
#Make this key argument
method = 3




# These can be input automatically from configurationn file
##################1701 Analysis##############

#raw_dir = '/Users/k3jackson/NBP1701/data/ctd_raw/'
#ssscc_file = '/Users/k3jackson/NBP1701/data/ssscc3.csv'
#time_dir = '/Users/k3jackson/NBP1701/data/time/'
#btl_dir = '/Users/k3jackson/NBP1701/data/bottle/'

##################1707 Analysis################
#raw_dir = '/Users/k3jackson/NPD1707/raw/'
#ssscc_file = '/Users/k3jackson/NPD1707/ssscc.csv'
#time_dir = '/Users/k3jackson/NPD1707/time/'
#btl_dir = '/Users/k3jackson/NPD1707/bottle/'

##################1706 Analysis################
raw_dir = '/Users/k3jackson/p06e/data/raw/'
ssscc_file = '/Users/k3jackson/p06e/data/ssscc.csv'
time_dir = '/Users/k3jackson/p06e/data/time/'
btl_dir = '/Users/k3jackson/p06e/data/bottle/'

#sal_col='CTDSAL'
#t_col='CTDTMP1'
#p_col='CTDPRS'
#lon_col='GPSLON'
#lat_col='GPSLAT'
#sal_btl_col='CTDSAL'
#t_btl_col='CTDTMP1'
#p_btl_col='CTDPRS'
#dov_col = 'CTDOXYVOLTS'
#lat_btl_col='GPSLAT'
#lon_btl_col='GPSLON'
#oxy_btl_col='CTDOXY1'
#dov_btl_col='CTDOXYVOLTS'


with open(ssscc_file, 'r') as filename:
    ssscc = [line.strip() for line in filename]

dataframe_concat = pd.DataFrame()

ssscc = ['02301']

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


    btl_data_write, btl_data_fit = oxy_fitting_ver2.oxy_fit(time_data,btl_data,stn_cst,hexfile,xmlfile)
    print('COMPLETED SSSCC:',stn_cst)
    dataframe_concat = pd.concat([dataframe_concat,btl_data_fit])
    
df=dataframe_concat
df['BTL_O'] = df['OXYGEN']-df['CTDOXY']
    
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
cm = ax.scatter(df['BTL_O'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
ax.set_xlim(-10,10)
ax.set_title('OXYGEN-CTDOXY vs CTDPRS')
ax.set_xlabel('CTDOXY Residual (umol/kg)')
ax.set_ylabel('Pressure (dbar)')
cbar = fig.colorbar(cm)
cbar.set_label('Station Number')