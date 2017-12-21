
# coding: utf-8

# # THIS IS ANALYSIS OF 1707 (PO6E) based upon NOAA's oxy fitting
# 
# ## Density Matching Instead of Interpolation

# In[1]:


import sys
sys.path.append('/Users/k3jackson/p06e/ctd_proc')
sys.path.append('/Users/k3jackson/p06e/Kennycode/')
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy
import numpy as np
import libODF_process_ctd as process_ctd
import libODF_sbe_reader as sbe_rd
import libODF_sbe_equations_dict as sbe_eq
import gsw
import pandas as pd



## In[2]:
#
#
#raw_dir = '/Users/k3jackson/NPD1707/raw/'
#ssscc_file = '/Users/k3jackson/NPD1707/ssscc.csv'
#time_dir = '/Users/k3jackson/NPD1707/time/'
#btl_dir = '/Users/k3jackson/NPD1707/bottle/'
#
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
#
#method = 3
#
#dataframe_concat = pd.DataFrame()
#ssscc = []
#with open(ssscc_file, 'r') as filename:
#    ssscc = [line.strip() for line in filename]
#for cast in range(len(ssscc)):
#    stn_nbr = int(ssscc[cast][0:3])
#    cst_nbr = int(ssscc[cast][-2:])
#    stn_cst = str(ssscc[cast])
#
#    hexfile = raw_dir+stn_cst+'.hex'
#    xmlfile = raw_dir+stn_cst+'.XMLCON'
#    
#    time_file = time_dir+stn_cst+'_time.pkl'
#    btl_file = btl_dir+stn_cst+'_btl_mean.pkl'
#    
#    time_data = process_ctd.dataToNDarray(time_file,float,True,',',1)
#    time_data = pd.DataFrame.from_records(time_data)  
#    
#    btl_data = process_ctd.dataToNDarray(btl_file,float,True,',',0) 
#    btl_data = pd.DataFrame.from_records(btl_data)
#

def oxy_fit(time_data,btl_data,ssscc,hexfile,xmlfile,method=3,sal_col='CTDSAL',\
            t_col='CTDTMP1',p_col='CTDPRS',lon_col='GPSLON',lat_col='GPSLAT',\
            dov_col='CTDOXYVOLTS',sal_btl_col='CTDSAL',t_btl_col='CTDTMP1',\
            p_btl_col='CTDPRS',lat_btl_col='GPSLAT',lon_btl_col='GPSLON',\
            dov_btl_col='CTDOXYVOLTS',oxy_btl_col='CTDOXY1'):
    
    ##NORMALLY INT#####
    stn_nbr = int(ssscc[0:3])
    cst_nbr = int(ssscc[-2:])
#    stn_cst = str(ssscc)
    
    time_data['CT_time'] = gsw.CT_from_t(time_data[sal_col],time_data[t_col],time_data[p_col])
    time_data['SA_time'] = gsw.SA_from_SP(time_data[sal_col],time_data[p_col],time_data[lon_col],time_data[lat_col])
    time_data['sigma0_time'] = gsw.sigma0(time_data['SA_time'],time_data['CT_time'])
    time_data_clean = process_ctd.binning_df(time_data)
    time_data_clean = time_data_clean.dropna(how='all')
    time_data_clean.index = pd.RangeIndex(len(time_data_clean.index))#reindex
    #bin dataset by pressure

    #pressure_seq_data = process_ctd.pressure_sequence(filename_base, p_col, timedate, 2.0, -1.0, 0.0, 'down', int(sample_rate), int(search_time), time_data)
    time_data_clean['CTDOXYVOLTS_time'] = time_data_clean[dov_col] #Added to remove confusion between btl and time data
    time_data_clean['TMP_KELVIN_CTD'] = time_data_clean[t_col] + 273.15
    time_data_clean['PRESSURE_CTD'] = time_data_clean[p_col]
    time_data_clean['TEMPERATURE_CTD'] = time_data_clean[t_col]
    time_data_clean['SALINITY_CTD'] = time_data_clean[sal_col]
    
    time_data_clean = time_data_clean.sort_values(by='sigma0_time')
    
    btl_data_clean = btl_data.dropna(how='all')
    btl_data_clean['BTLNBR'] = btl_data['btl_fire_num'] #May not be best way to get bottle number
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index)) #reindex
    btl_data_clean['CT_btl'] = gsw.CT_from_t(btl_data_clean[sal_btl_col],btl_data_clean[t_btl_col],btl_data_clean[p_btl_col])
    btl_data_clean['SA_btl'] = gsw.SA_from_SP(btl_data_clean[sal_btl_col],btl_data_clean[p_btl_col],btl_data_clean[lon_btl_col],btl_data_clean[lat_btl_col])
    btl_data_clean['sigma0_btl'] = gsw.sigma0(btl_data_clean['SA_btl'],btl_data_clean['CT_btl'])

    #Sort bottles for matching
    btl_data_clean = btl_data_clean.sort_values(by='sigma0_btl') #reindex
    
    doxyv_btl = np.diff(btl_data_clean['CTDOXYVOLTS'])
    dt_btl = np.diff(btl_data_clean['scan_datetime'])
    dv_dt_btl = doxyv_btl/dt_btl
    dv_dt_conv_btl = np.convolve(dv_dt_btl,[0.5,0.5])
    btl_data_clean['dv_dt_conv_btl'] = dv_dt_conv_btl

    
    #plt.plot(dv_dt_conv_btl)
    
    #calculate dV/dT and filter signal
    
    doxyv = np.diff(time_data_clean['CTDOXYVOLTS_time'])
    dt = np.diff(time_data_clean['scan_datetime'])
    
    dv_dt = doxyv/dt
    dv_dt_conv = np.convolve(dv_dt,[0.5,0.5])
    
    a = 1
    windowsize = 5
    b = (1/windowsize)*np.ones(windowsize)
    
    #filt = scipy.signal.lfilter(b,a,dv_dt_conv)
    filtfilt = scipy.signal.filtfilt(b,a,dv_dt_conv)
    
    #plt.plot(dv_dt_conv,'r',label='original convolved signal')
    #plt.plot(filt,label='filtered 1 direction')
    #plt.plot(filtfilt,'g',label='filtered 2 directions')
    #plt.legend()
    #plt.title('dv_dt filtering comparison')
    
    time_data_clean['dv_dt_time'] = filtfilt
    
    all_sigma0 = time_data_clean['sigma0_time'].append(btl_data_clean['sigma0_btl'])
    fsg = pd.Series(np.min(all_sigma0)-1e-4)
    lsg = pd.Series(np.max(all_sigma0)+1e-4)
    new_sigma0 = fsg.append(time_data_clean['sigma0_time'])
    new_sigma0 = new_sigma0.append(lsg)
    np.size(new_sigma0)

    new_pressure = pd.Series(time_data_clean['CTDPRS'].iloc[0])
    new_pressure = new_pressure.append(time_data_clean['CTDPRS'])
    new_pressure = new_pressure.append(pd.Series(time_data_clean['CTDPRS'].iloc[-1]))
    len(new_pressure)
    
    new_sal = pd.Series(time_data_clean['CTDSAL'].iloc[0])
    new_sal = new_sal.append(time_data_clean['CTDSAL'])
    new_sal = new_sal.append(pd.Series(time_data_clean['CTDSAL'].iloc[-1]))
    np.size(new_sal)

    new_temp = pd.Series(time_data_clean['CTDTMP1'].iloc[0])
    new_temp = new_temp.append(time_data_clean['CTDTMP1'])
    new_temp = new_temp.append(pd.Series(time_data_clean['CTDTMP1'].iloc[-1]))
    np.size(new_temp)

    new_oxyvo = pd.Series(time_data_clean['CTDOXYVOLTS_time'].iloc[0])
    new_oxyvo = new_oxyvo.append(time_data_clean['CTDOXYVOLTS_time'])
    new_oxyvo = new_oxyvo.append(pd.Series(time_data_clean['CTDOXYVOLTS_time'].iloc[-1]))
    np.size(new_oxyvo)

    x = np.arange(np.size(time_data_clean['sigma0_time']))
    x_inter = np.arange(np.size(new_sigma0))
    inter_sigma2 = scipy.interpolate.interp1d(x_inter,new_sigma0)
    new_x = np.linspace(0,np.max(x_inter),np.size(time_data_clean['sigma0_time'])) #Make linspace of len of time_data_cleanto get new sigmas and add to dataframe
    new_sigma0=inter_sigma2(new_x)

    inter_pres2 = scipy.interpolate.interp1d(x_inter,new_pressure)
    inter_sal2 = scipy.interpolate.interp1d(x_inter,new_sal)
    inter_temp2 = scipy.interpolate.interp1d(x_inter,new_temp)
    inter_oxyvo2 = scipy.interpolate.interp1d(x_inter,new_oxyvo)
    
    #plt.plot(new_sigma0)
    
    dpres = inter_pres2(new_x)
    dsal = inter_sal2(new_x)
    dtemp = inter_temp2(new_x)
    doxyvo = inter_oxyvo2(new_x)
    
    #dv/dt Interpolation
    btl_len_x = np.arange(np.size(btl_data_clean['dv_dt_conv_btl']))
    dv_dt_inter_x = np.arange(np.size(time_data_clean['dv_dt_time']))
    dv_dt_inter = np.interp(btl_len_x,dv_dt_inter_x,time_data_clean['dv_dt_time'])
    
    time_data_clean['sigma0_time'] = new_sigma0
    time_data_clean['PRESSURE_CTD'] = dpres
    time_data_clean['SALINITY_CTD'] = dsal
    time_data_clean['TEMPERATURE_CTD'] = dtemp
    time_data_clean['CTDOXYVOLTS_time'] = doxyvo

    btl_data_clean = btl_data.dropna(how='all')
    btl_data_clean['BTLNBR'] = btl_data['btl_fire_num'] #May not be best way to get bottle number
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index)) #reindex
    btl_data_clean['CT_btl'] = gsw.CT_from_t(btl_data_clean[sal_btl_col],btl_data_clean[t_btl_col],btl_data_clean[p_btl_col])
    btl_data_clean['SA_btl'] = gsw.SA_from_SP(btl_data_clean[sal_btl_col],btl_data_clean[p_btl_col],btl_data_clean[lon_btl_col],btl_data_clean[lat_btl_col])
    btl_data_clean['sigma0_btl'] = gsw.sigma0(btl_data_clean['SA_btl'],btl_data_clean['CT_btl'])

    #Sort bottles for matching
    btl_data_clean = btl_data_clean.sort_values(by='sigma0_btl') #reindex

    matched_df = pd.merge_asof(btl_data_clean,time_data_clean,left_on='sigma0_btl',right_on='sigma0_time',direction='nearest')
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index))

    time_data_matched = matched_df
    
    time_data_matched['STNNBR']=stn_nbr
    time_data_matched['CASTNO']=cst_nbr

    btl_data_clean['STNNBR']=stn_nbr
    btl_data_clean['CASTNO']=cst_nbr
    
    time_data_matched['dv_dt_time'] = dv_dt_inter*-1
    
    btl_data_clean['OS'] = sbe_eq.OxSol(btl_data_clean['CTDTMP1'],btl_data_clean['CTDSAL'])
    time_data_matched['OS'] = sbe_eq.OxSol(time_data_matched['TEMPERATURE_CTD'],time_data_matched['SALINITY_CTD'])
    
    coef0 = None
    
    # GET Initial Guess from hex file
    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1, 'ranking': 5, 'column': 'CTDOXYVOLTS', 'sensor_info': rawConfig['Sensors'][i]}

    if coef0 is None:
        coef0 = [oxy_meta['sensor_info']['Soc'], oxy_meta['sensor_info']['offset'], oxy_meta['sensor_info']['A'], oxy_meta['sensor_info']['B'], oxy_meta['sensor_info']['C'], oxy_meta['sensor_info']['E']]
        print(coef0)

    Tau20=oxy_meta['sensor_info']['Tau20']
    Tcorr = oxy_meta['sensor_info']['Tcor']
    
    coef0 =[oxy_meta['sensor_info']['Soc'],oxy_meta['sensor_info']['offset'],oxy_meta['sensor_info']['Tau20'],oxy_meta['sensor_info']['Tcor'],oxy_meta['sensor_info']['E']]
    
    # D1 and D2 are fixed here

    cc=[1.92634e-4,-4.64803e-2]

    #cc(1)=1.92634e-4;
    #cc(2)=-4.64803e-2;
    
    btl_data_clean['dv_dt_conv_btl'] = dv_dt_conv_btl
    
    btl_data_clean['NOAA_oxy_mlL'] =coef0[0]*(btl_data_clean['CTDOXYVOLTS']+coef0[1]+Tau20*np.exp(cc[0]*btl_data_clean['CTDPRS']+cc[0]*btl_data_clean['CTDTMP1'])
             *btl_data_clean['dv_dt_conv_btl'])*btl_data_clean['OS']*np.exp(Tcorr*btl_data_clean['CTDTMP1'])*np.exp((coef0[4]*btl_data_clean['CTDPRS'])/(btl_data_clean['CTDTMP1']+273.15))

    time_data_matched['NOAA_oxy_mlL'] =coef0[0]*(time_data_matched['CTDOXYVOLTS_y']+coef0[1]+Tau20*np.exp(cc[0]*time_data_matched['CTDPRS_y']+cc[0]*time_data_matched['TEMPERATURE_CTD'])
             *time_data_matched['dv_dt_time'])*time_data_matched['OS']*np.exp(Tcorr*time_data_matched['TEMPERATURE_CTD'])*np.exp((coef0[4]*time_data_matched['CTDPRS_y'])/(time_data_matched['TEMPERATURE_CTD']+273.15))
    
    eps=1e-5
    wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
    wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]
    wgt = scipy.interpolate.interp1d(wrow1,wrow2)
    
    time_data_matched['weights'] = wgt(time_data_matched['CTDPRS_y'])
    
    coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef0,args=(time_data_matched,btl_data_clean,method))
    
    new_OXY = coef[0]*(time_data_matched['CTDOXYVOLTS_y']+coef[1]+coef[2]*np.exp(cc[0]*time_data_matched['CTDPRS_y']+cc[0]*time_data_matched['TEMPERATURE_CTD'])    *time_data_matched['dv_dt_time'])*time_data_matched['OS']*np.exp(coef[3]*time_data_matched['TEMPERATURE_CTD'])    *np.exp((coef[4]*time_data_matched['CTDPRS_y'])/(time_data_matched['TEMPERATURE_CTD']+273.15))
    
    time_data_matched['Fitted_OXY_uMOLKG_time'] = new_OXY * 44660 / (time_data_matched['sigma0_time'] + 1000)
    btl_data_clean['OXY_uMOLKG_btl'] = btl_data_clean[oxy_btl_col] * 44660 / (btl_data_clean['sigma0_btl'] + 1000)
    btl_data_clean['CTDOXY'] = time_data_matched['Fitted_OXY_uMOLKG_time']
    btl_data_output = btl_data_clean
    
    btl_data_output['residual'] = np.abs(btl_data_clean['OXY_uMOLKG_btl']-time_data_matched['Fitted_OXY_uMOLKG_time'])
    cutoff = np.std(btl_data_output['residual'])*2.8
    btl_data_output['index'][btl_data_output['residual']<=cutoff]
    btl_data_clean=btl_data_output[btl_data_output['residual']<=cutoff]
    thrown_values = btl_data_output[btl_data_output['residual']>cutoff]
    bad_values = btl_data_output[btl_data_output['residual']>cutoff]
    btl_data_clean
    
        #Repeat fitting until there are no outliers
    while not thrown_values.empty:
        matched_df = pd.merge_asof(btl_data_clean,time_data_clean,left_on='sigma0_btl',right_on='sigma0_time',direction='nearest')
        btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index)) #reindex DF
        time_data_matched = matched_df
        #RE-interpolate dv_dt
        btl_len_x = np.arange(np.size(btl_data_clean['dv_dt_conv_btl']))
        dv_dt_inter_x = np.arange(np.size(time_data_clean['dv_dt_time']))
        dv_dt_inter = np.interp(btl_len_x,dv_dt_inter_x,time_data_clean['dv_dt_time'])
        time_data_matched['dv_dt_time'] = dv_dt_inter*-1
        #Calculate SOL
        
        btl_data_clean['OS'] = sbe_eq.OxSol(btl_data_clean['CTDTMP1'],btl_data_clean['CTDSAL'])
        time_data_matched['OS'] = sbe_eq.OxSol(time_data_matched['TEMPERATURE_CTD'],time_data_matched['SALINITY_CTD'])
        
        #Recalculate new oxy and coefficients using new values and previously determined coef as initial guess
        
        btl_data_clean['NOAA_oxy_mlL'] =coef[0]*(btl_data_clean['CTDOXYVOLTS']+coef[1]+Tau20*np.exp(cc[0]*btl_data_clean['CTDPRS']+cc[0]*btl_data_clean['CTDTMP1'])
        *btl_data_clean['dv_dt_conv_btl'])*btl_data_clean['OS']*np.exp(Tcorr*btl_data_clean['CTDTMP1'])*np.exp((coef[4]*btl_data_clean['CTDPRS'])/(btl_data_clean['CTDTMP1']+273.15))
        
        time_data_matched['NOAA_oxy_mlL'] =coef[0]*(time_data_matched['CTDOXYVOLTS_y']+coef[1]+Tau20*np.exp(cc[0]*time_data_matched['CTDPRS_y']+cc[0]*time_data_matched['TEMPERATURE_CTD'])
        *time_data_matched['dv_dt_time'])*time_data_matched['OS']*np.exp(Tcorr*time_data_matched['TEMPERATURE_CTD'])*np.exp((coef[4]*time_data_matched['CTDPRS_y'])/(time_data_matched['TEMPERATURE_CTD']+273.15))
        
        coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef,args=(time_data_matched,btl_data_clean,method))
        #coef = fit_ctd.find_oxy_coef(btl_data_clean[oxy_btl_col], \
        #btl_data_clean[p_btl_col],btl_data_clean[t_btl_col], \
        #btl_data_clean[sal_btl_col],btl_data_clean[dov_btl_col],hexfile,xmlfile,coef)

        #Apply new coef to time_data
        
        new_OXY = coef[0]*(time_data_matched['CTDOXYVOLTS_y']+coef[1]+coef[2]*np.exp(cc[0]*time_data_matched['CTDPRS_y']+cc[0]*time_data_matched['TEMPERATURE_CTD'])*time_data_matched['dv_dt_time'])*time_data_matched['OS']*np.exp(coef[3]*time_data_matched['TEMPERATURE_CTD'])*np.exp((coef[4]*time_data_matched['CTDPRS_y'])/(time_data_matched['TEMPERATURE_CTD']+273.15))
        #time_data_matched['Fitted_OXY_mLL'] = fit_ctd.oxy_dict(coef,\
        #time_data_matched['PRESSURE_CTD'], time_data_matched['TMP_KELVIN_CTD'], \
        #time_data_matched['TEMPERATURE_CTD'], time_data_matched['SALINITY_CTD'], \
        #time_data_matched['CTDOXYVOLTS_time'])

        #Calculate Oxy in uMol/Kg
        time_data_matched['Fitted_OXY_uMOLKG_time'] = new_OXY * 44660 / (time_data_matched['sigma0_time'] + 1000)
        btl_data_clean['OXY_uMOLKG_btl'] = btl_data_clean[oxy_btl_col] * 44660 / (btl_data_clean['sigma0_btl'] + 1000)

        plt.plot(btl_data_clean['OXY_uMOLKG_btl'],btl_data_clean['CTDPRS'],label='btl data')
        plt.plot(time_data_matched['Fitted_OXY_uMOLKG_time'],time_data_matched['CTDPRS_y'],'g',label='time data')
        plt.legend()
        plt.xlabel('Oxygen (uMol/KG)')
        plt.ylabel('Pressure (dBar)')
        plt.show()

        btl_data_clean['OXYGEN'] = btl_data_clean['OXY_uMOLKG_btl']
        btl_data_clean['CTDOXY'] = time_data_matched['Fitted_OXY_uMOLKG_time']

        btl_data_output = btl_data_clean
        btl_data_clean['residual'] = np.abs(btl_data_output['OXYGEN']-btl_data_output['CTDOXY'])
        std_res = np.std(btl_data_output['residual'])
        cutoff = np.std(btl_data_output['residual'])*2.8
        
        btl_data_clean=btl_data_clean[btl_data_output['residual']<=cutoff]
        bad_values2 = btl_data_output[btl_data_output['residual']>cutoff]
        thrown_values = btl_data_clean[btl_data_output['residual']>cutoff]
        bad_values = pd.concat([bad_values,bad_values2])
        print('NEW STANDARD DEVIATION IS:',std_res)
    #Add Station and Cast Numbers to DataFrame
    time_data_matched['STNNBR']=stn_nbr
    time_data_matched['CASTNO']=cst_nbr

    btl_data_clean['STNNBR']=stn_nbr
    btl_data_clean['CASTNO']=cst_nbr
    
    plt.plot(btl_data_clean['residual'],btl_data_clean['CTDPRS']*-1,'bx')
    plt.xlim(xmax=10)
    plt.show()
    
    bad_values['CTDOXY_FLAG_W'] = 4
    btl_data_clean['CTDOXY_FLAG_W'] = 2
    btl_data_clean = pd.concat([btl_data_clean,bad_values])
    btl_data_clean=btl_data_clean.sort_values(by='BTLNBR')
    
    btl_data_clean = btl_data_clean.reset_index(drop=True)
    btl_data_write = pd.DataFrame()
    btl_data_write['STNNBR'] = btl_data_clean['STNNBR'].values
    btl_data_write['CASTNO'] = btl_data_clean['CASTNO'].values
    btl_data_write['SAMPNO'] = btl_data_clean['BTLNBR'].values
    btl_data_write['CTDOXY'] = btl_data_clean['CTDOXY'].values
    btl_data_write['CTDOXY_FLAG_W'] = btl_data_clean['CTDOXY_FLAG_W'].values
    
    return btl_data_write, btl_data_clean

def oxygen_cal_ml(coef0,time_data,btl_data,switch):

    """"

    NOAA's oxygen fitting routine using the equation:
    OXY(ml/L) = SOC * (doxy_volts+Voffset+Tau20*exp(cc1*PRES+cc2*TEMP)*dvdt) \
                *os*exp(Tcorr*TEMP)*exp(E*PRESS/TEMP_K)

    coef0s:
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E

    cc[0] = D1
    cc[1] = D2
    
    """

    cc = [1.92634e-4,-4.64803e-2]
    #ORIGINAL CODE
    #time_data['NOAA_oxy_mlL'] =coef0[0]*(time_data['CTDOXYVOLTS_y']+\
    #coef0[1]+coef0[2]*np.exp(cc[0]*time_data['CTDPRS_y']+cc[0]*time_data['TEMPERATURE_CTD'])\
    #*time_data['dv_dt_time'])*time_data['OS']*np.exp(coef0[3]*time_data['TEMPERATURE_CTD'])\
    #*np.exp((coef0[4]*time_data['CTDPRS_y'])/(time_data['TEMPERATURE_CTD']+273.15))

    #MODIFIED CODE
    time_data['NOAA_oxy_mlL'] =coef0[0]*(time_data['CTDOXYVOLTS_time']+\
    coef0[1]+coef0[2]*np.exp(cc[0]*time_data['PRESSURE_CTD']+cc[0]*time_data['TEMPERATURE_CTD'])\
    *time_data['dv_dt_time'])*time_data['OS']*np.exp(coef0[3]*time_data['TEMPERATURE_CTD'])\
    *np.exp((coef0[4]*time_data['PRESSURE_CTD'])/(time_data['TEMPERATURE_CTD']+273.15))

    #Weight Determination
    if switch == 1:
        eps=1e-5
        wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
        wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]
        wgt = scipy.interpolate.interp1d(wrow1,wrow2)

        #Original Code
        #time_data['weights'] = wgt(time_data['CTDPRS_y'])

        #Modified CODE
        time_data['weights'] = wgt(time_data['PRESSURE_CTD'])

        #resid = np.sum((time_data['weights']*(btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL']))/(np.sum(time_data['weights'])**2))
        #resid = (time_data['weights']*(btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL']))/(np.sum(time_data['weights'])**2) Working
        resid = (time_data['weights']*(btl_data['NOAA_oxy_mlL']-time_data['NOAA_oxy_mlL']))/(np.sum(time_data['weights'])**2) 
    elif switch ==2:
        #L2 Norm
        resid = (btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL'])**2
    elif switch ==3:
        #ODF residuals
        resid = np.sqrt(((btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL'])**2)/(np.std(time_data['NOAA_oxy_mlL'])**2))

    #elif switch ==4:

    return resid


def match_series_len(df,col_name):
    ''' This code simply makes a copy of the first and last elements of a series and appends them
        to the beginning and end of a series respectively

    Inputs
        df: DataFrame
        col_name: Name of column that needs to be worked

    Outputs
        new_series: appended data series

    '''

    new_series = pd.Series(df[col_name].iloc[0])
    new_series = new_series.append(df['CTDPRS'])
    new_series = new_series.append(pd.Series(df[col_name].iloc[-1]))

    return new_series