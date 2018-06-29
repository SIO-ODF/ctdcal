#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:04:19 2018

@author: k3jackson
"""


import scipy
import numpy as np
import sys
sys.path.append('ctdcal/')
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import ctdcal.sbe_reader as sbe_rd
import ctdcal.sbe_equations_dict as sbe_eq
import gsw
import pandas as pd
import csv

#Line 342 module isopycnals

def oxy_loader(oxyfile):
    
    f = open(oxyfile,newline='')
    oxyF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
    oxy_Array = []
    
    for row in oxyF:
        if len(row) > 9:
            row = row[:9]
        oxy_Array.append(row)
    
    params = oxy_Array[0]
    
    del oxy_Array[0]
        
    header = ['STNNO','CASTNO','BOTTLENO_OXY','FLASKNO','TITR_VOL','TITR_TEMP','DRAW_TEMP','TITER_TIME','END_VOLTS']
    
    O2 = np.array(oxy_Array)
    df = pd.DataFrame(O2,columns=header)
    
    f.close()
    
    df = df.apply(pd.to_numeric, errors='ignore')
    
    # Remove "Dummy Data"
    #df['BOTTLENO_OXY'] = df['BOTTLENO_OXY'].astype(int,errors='ignore')
    df = df[df['BOTTLENO_OXY']!=99]
    
    # Remove "ABORTED DATA"
    
    df = df[df['TITR_VOL'] > 0]
    
    # Sort and reindex
    df = df.sort_values('BOTTLENO_OXY')
    #df = df.reset_index()
    
    # Get necessary columns for output
    df['STNNO'] = df['STNNO'].astype(str)
    df['CASTNO'] = df['CASTNO'].astype(str)
    df['FLASKNO']= df['FLASKNO'].astype(str)
    
    df['SSSCC_OXY'] = df['STNNO']+'0'+df['CASTNO']
    
    # Reindex dataframe to have 36 places
    DF = pd.DataFrame(data=np.arange(1,37),columns=['BOTTLENO_OXY'],index=range(1,37))
    DF = DF.merge(df,on="BOTTLENO_OXY",how='outer')
    DF = DF.set_index(np.arange(1,37))

    
    return DF ,params

def flask_load(flaskfile='data/oxygen/o2flasks.vol',skip_rows=12):
    """Load information from flask.vol file
    
    """
## Convert to Dataframe
    with open(flaskfile, 'r') as f:
        flasks = {}
        for l in f:
            if 'Volume' not in l:
                if l.strip().startswith("#"):
                    continue
                row = l.strip().split()
                flasks[str(row[0])] = float(row[1])
                
    flasks = pd.DataFrame.from_dict(flasks,orient='index')
    flasks = flasks.rename(columns={0:'FLASK_VOL'})
    flasks['FLASKNO'] = flasks.index.values
    flasks = flasks[['FLASKNO','FLASK_VOL']]
    return flasks
    
#    f = open(flaskfile,newline='')
#    flaskF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
#
#    flask_array = []
#    for row in flaskF:
#        flask_array.append(row)
#    
#    f.close
#    del flask_array[0:skip_rows]
#    
#    df = pd.DataFrame(flask_array)
#    df['FLASKNO'] = df[0] 
#    df['FLASK_VOL'] = df[1]
#    df = df[['FLASKNO','FLASK_VOL']]
#    df['FLASKNO'] = df['FLASKNO'].astype(float,errors='ignore')
    
#    return df

def match_flask(df,flask_df,ssscc_col='SSSCC',vol_col='FLASK_VOL',num_col='FLASKNO',t=20, glass="borosilicate"):
    ### THIS CODE MAY NEED AN ADAPTATION FOR A NON DF INPUT. E.G. PANDAS SERIES
    
    ### NOTE: LOOK AT ANDREW DICKSON'S SOP 13: GRAVIMETRIC CALIBRATION OF A VOLUME CONTAINED USING WATER!!!
    ### These formulas/values may be out of date
    _coef = {
            "borosilicate": 0.00001,
            "soft": 0.000025,
            }
    _t = 20
    coef = _coef[glass]
    merge_df = df.merge(flask_df,how='outer')
    merge_df = merge_df[merge_df[ssscc_col].notnull()]
    #merge_df = merge_df.reset_index()
#    
#    if isinstance(t,int):
#        merge_df[vol_col] = merge_df[vol_col] * (1.0 + coef * (t - _t))
#        
#    elif isinstance(t,pd.core.series.Series):
#        t.index = merge_df.index.values
#        merge_df[vol_col] = merge_df[vol_col] * (1.0 + coef * (t - _t))
    ### CHANGE TO A MERGE/JOIN FUNCTION
#    flask_array = []
#    for x in flask:
#        if x != np.nan:
#            #flask_vol = float(flask_df[vol_col][flask_df[num_col]==x].values[0])
#            flask_vol = flask_df[vol_col][flask_df[num_col]==x]
#            flask_vol = float(flask_vol[0].values[0])
#        
#            coef = _coef[glass]
#            flask_vol = flask_vol * (1.0 + coef * (t - _t));
#            flask_array.append(flask_vol)
#        elif x ==np.nan:
#            flask_array.append(np.nan)
    merge_df = merge_df.sort_values(by='master_index')
    merge_df = merge_df.set_index(df.index.values)    
    
    return merge_df

def gather_oxy_params(oxyfile):

    f = open(oxyfile,newline='')
    oxyF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
    oxy_Array = []
    
    for row in oxyF:
        oxy_Array.append(row)
    
    f.close()
    
    params = oxy_Array[0]   
    df = pd.DataFrame(params)
    df = df.transpose() 
    
    df = df[[0,1,2,3,4,5]]
    
# If column names are needed later on for the parameters DF    
#    cols = ['TITR','BLANK','KIO3_N','KIO3_V','KIO3_T','THIO_T']
#    
#    df.rename(columns={0:cols[0],1:cols[1],2:cols[2],3:cols[3],4:cols[4],5:cols[5]})
    
    return df

def gather_all_oxy_params(ssscc,oxyfile_prefix='data/oxygen/',oxyfile_postfix=''):
    
    df_data_all = pd.DataFrame()
    
    for x in ssscc:
        oxyfile = oxyfile_prefix + x + oxyfile_postfix
        params = gather_oxy_params(oxyfile)
        params['SSSCC'] = x
        
        df_data_all = pd.concat([df_data_all,params])
        
    df_data_all = df_data_all.reset_index()
    df_data_all = df_data_all.set_index('SSSCC')
    return df_data_all

def thio_n_calc(df,params,ssscc_col='SSSCC',thio_col = 'thio_n'):#(titr, blank, kio3_n, kio3_v, kio3_t, thio_t):

    params = params.loc[df[ssscc_col]]
    params = params.apply(pd.to_numeric)
    
    titr   = params[0]
    blank  = params[1]
    kio3_n = params[2]
    kio3_v = params[3]
    kio3_t = params[4]
    thio_t = params[5]
    
    rho_stp  = fit_ctd.rho_t(20)
    rho_kio  = fit_ctd.rho_t(kio3_t)
    rho_thio = fit_ctd.rho_t(thio_t)

    kio3_v_20c = kio3_v * (rho_kio/rho_stp)
    thio_v_20c = titr * (rho_thio/rho_stp) - blank

    thio_n = kio3_v_20c * kio3_n / thio_v_20c
    thio_n = pd.Series(thio_n.values,index=df.index.values)
    df[thio_col] = thio_n
    return df

def titr_20_calc(titr, titr_temp):
    rho_titr = fit_ctd.rho_t(titr_temp)
    rho_stp = fit_ctd.rho_t(20)

    #convert the titr ammount to 20c equivalent
    titr_20c = titr * rho_titr/rho_stp
    return titr_20c

def calculate_bottle_oxygen(df,ssscc,ssscc_col='SSSCC',vol_col='FLASK_VOL',titr_temp_col='TITR_TEMP',titr_col='titr_20C',thio_n_col='thio_n',d_temp='DRAW_TEMP'):
    

    
    params = gather_all_oxy_params(ssscc)
    
    df = thio_n_calc(df,params)
    
    df[titr_col] = titr_20_calc(df['TITR_VOL'],df['TITR_TEMP'])
    
    flask_df = flask_load()
    
    df = match_flask(df,flask_df,t=df[d_temp])
    
    params = params.loc[df[ssscc_col]]
    params = params.apply(pd.to_numeric)    
    
    blank  = params[1]
    
    blank = blank.reset_index()
    blank = blank.set_index(df.index.values)
    blank = blank[1]

    oxy_mlL = oxygen_eq(df[titr_col],blank,df[thio_n_col],df[vol_col])
#    oxy_mlL = (((titr_20c - blank) * thio_n * 5.598 - 0.0017)/((df[vol_col] - 2.0) * 0.001))
    
    
    return oxy_mlL

def oxygen_eq(titr,blank,thio_n,flask_vol):
    
    oxy_mlL = (((titr - blank) * thio_n * 5.598 - 0.0017)/((flask_vol - 2.0) * 0.001))
    
    return oxy_mlL

def calibrate_oxygen(df_time,df_btl,ssscc,method=3,raw_dir='data/raw/'):
    
    btl_concat = pd.DataFrame()
    time_concat = pd.DataFrame()
    for x in ssscc:

        hexfile = raw_dir+x+'.hex'
        xmlfile = raw_dir+x+'.XMLCON'
        
        time_data = df_time[df_time['SSSCC'] == x]
        btl_data = df_btl[df_btl['SSSCC'] == x]
        coef = oxy_fit(time_data,btl_data,x,hexfile,xmlfile,method=method)
        #btl_data_write, btl_data_fit = oxy_fitting.oxy_fit(time_data,btl_data,x,hexfile,xmlfile,method=method)
        print('COMPLETED SSSCC:',x)
        #dataframe_concat = pd.concat([dataframe_concat,btl_data_fit])
        
        btl_data['OXYGEN'] = apply_oxy_coef(btl_data,coef)
        time_data['CTDOXY'] = apply_oxy_coef(time_data,coef)
    
        btl_concat = pd.concat([btl_concat,btl_data])
        time_concat = pd.concat([time_concat,time_data])
          
    return btl_concat, time_concat

def oxy_from_ctd_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc=[1.92634e-4,-4.64803e-2]):
    """
    
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E
    """
    
    
   # ox=c(1)*(doxvo+c(2)+c(3)*exp(cc(1)*dpres+cc(2)*dtepr).*dsdvdt).*os.*exp(c(4)*dtepr).*exp(c(5)*dpres./(273.15+dtepr));
    oxy_mlL = coef[0] * (oxyvolts + coef[1] + coef[2] * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt) * os \
            * np.exp(coef[3] * temp) \
            * np.exp((coef[4] * pressure) \
            / (temp+273.15))
           
    
    return oxy_mlL

def SB_oxy_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc=[1.92634e-4,-4.64803e-2]):
    """    
    Oxygen (ml/l) = Soc * (V + Voffset) * (1.0 + A * T + B * T + C * T ) * OxSol(T,S) * exp(E * P / K)
    Original equation from calib sheet dated 2014:
    
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = A
    coef[3] = B
    coef[4] = C
    coef[5] = E
    coef[6] = Tau20
        
    """    
    
    oxygen = (coef[0]) * ((oxyvolts + coef[1] + (coef[6] * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt))
                 * (1.0 + coef[2] * temp + coef[3] * np.power(temp,2) + coef[4] * np.power(temp,3) )
                 * os
                 * np.exp((coef[5] * pressure) / (temp+273.15)))
    
    return oxygen

def get_SB_coef(hexfile,xmlfile):
    
    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1, 
                        'ranking': 5, 'column': 'CTDOXYVOLTS', 
                        'sensor_info': rawConfig['Sensors'][i]}

     
    Soc = oxy_meta['sensor_info']['Soc'] # Oxygen slope
    Voff = oxy_meta['sensor_info']['offset'] # Sensor output offset voltage
    Tau20 = oxy_meta['sensor_info']['Tau20'] #Sensor time constant at 20 deg C and 1 Atm
    Tcorr = oxy_meta['sensor_info']['Tcor'] # Temperature correction
    E  = oxy_meta['sensor_info']['E'] #Compensation Coef for pressure effect on membrane permeability (Atkinson et al)
    A = oxy_meta['sensor_info']['A'] # Compensation Coef for temp effect on mem perm
    B = oxy_meta['sensor_info']['B'] # Compensation Coef for temp effect on mem perm
    C = oxy_meta['sensor_info']['C'] # Compensation Coef for temp effect on mem perm
    D = [oxy_meta['sensor_info']['D0'],oxy_meta['sensor_info']['D1'],oxy_meta['sensor_info']['D2']] # Press effect on time constant Usually not fitted for.
    
    #coef = {'Soc':Soc,'Voff':Voff,'Tau20':Tau20,'Tcorr':Tcorr,'E':E,'A':A,'B':B,'C':C,'D':D}
    
    # Make into an array in order to do keast squares fitting easier
#        coef0s:
#    coef[0] = Soc
#    coef[1] = Voffset
#    coef[2] = Tau20
#    coef[3] = Tcorr
#    coef[4] = E
#
#    cc[0] = D1
#    cc[1] = D2
    coef = [Soc,Voff,A,B,C,E,Tau20]
    
    return coef


def apply_oxy_coef(df,coef,cc=[1.92634e-4,-4.64803e-2],oxyvo_col='CTDOXYVOLTS',p_col='CTDPRS',t_col='CTDTMP1',sal_col='CTDSAL'):
    
    sigma0 = sigma0_from_CTD(df,sal_col=sal_col,t_col=t_col,p_col=p_col)
    
    df['OS'] = sbe_eq.OxSol(df[t_col],df[sal_col])
    df['dv_dt'] = calculate_dVdT(df)
    
    new_OXY = coef[0] * (df[oxyvo_col] + coef[1] + coef[2] \
            * np.exp(cc[0] * df[p_col] + cc[0] \
            * df[t_col]) \
            * df['dv_dt']) * df['OS'] \
            * np.exp(coef[3] * df[t_col]) \
            * np.exp((coef[4] * df[p_col]) \
            / (df[t_col]+273.15))
            
    oxy_uMolkg = oxy_ml_to_umolkg(new_OXY,sigma0)
    
    return oxy_uMolkg
    
def sigma_from_CTD(df,sal_col='CTDSAL',t_col='CTDTMP1',p_col='CTDPRS',lon_col='GPSLON',lat_col='GPSLAT',ref=0):
    
    CT = gsw.CT_from_t(df[sal_col],df[t_col],df[p_col])
    
    SA = gsw.SA_from_SP(df[sal_col],df[p_col],df[lon_col],df[lat_col])
    
    # Reference pressure in ref*1000 dBars
    if ref == 0:
        sigma = gsw.sigma0(SA,CT)
    elif ref == 1:
        sigma = gsw.sigma1(SA,CT)
    elif ref == 2:
        sigma = gsw.sigma2(SA,CT)
    elif ref == 3:
        sigma = gsw.sigma3(SA,CT)
    elif ref == 4:
        sigma = gsw.sigma4(SA,CT)
    
    return sigma

def os_umol_kg(SA,temp,press):
    """ 
    Oxygen solubility in umol/kg from Gordon and Garcia 1992
    """
    
    
    
    a0=5.80871
    a1=3.20291
    a2=4.17887
    a3=5.10006
    a4=-9.86643e-2
    a5=3.80369
    b0=-7.01577e-3
    b1=-7.70028e-3
    b2=-1.13864e-2
    b3=-9.51519e-3
    c0=-2.75915e-7
    
    
    # Calculate potential temp
    PT = gsw.pt0_from_t(SA,temp,press)
    
    TS= np.log((298.15-PT)/(273.15+PT))
    
    x= a0 + a1 * TS + a2 * TS**2 + a3 * TS**3 + a4 * TS**4 + a5 * TS**5 + SA * (b0 + b1 * TS + b2 * TS**2 + b3 * TS**3) + c0 * SA**2
    
    os_umol = np.exp(x)
    
    return os_umol

def calculate_dVdT(df,oxyvo_col='CTDOXYVOLTS',time_col='scan_datetime'):
#    
     doxyv = df[oxyvo_col].diff()
     dt = df[time_col].diff()
     dv_dt = doxyv/dt
     dv_dt_mean = np.mean(dv_dt[(dv_dt != np.inf) & (dv_dt !=-np.inf)])
     dv_dt = dv_dt.replace([np.inf,-np.inf],np.NaN)
     dv_dt = dv_dt.fillna(dv_dt_mean)
     
     #dv_dt_conv = np.convolve(dv_dt,[0.5,0.5])
     
     a = 1
     windowsize = 5
     b = (1/windowsize)*np.ones(windowsize)
     #filtered_dvdt = scipy.signal.filtfilt(b,a,dv_dt_conv)
     filtered_dvdt = scipy.signal.filtfilt(b,a,dv_dt)
     
     
     return filtered_dvdt
 
def merge_parameters(btl_df,time_df,l_param='sigma0_btl',r_param='sigma0_ctd'):
        
    merge_df = pd.merge_asof(btl_df,time_df,left_on=l_param,right_on=r_param,direction='nearest')
    
    # Clean up merged dataframe
    # These column names can be taken as a variable from the ini file
    keep_columns = ['CTDTMP1_y','CTDTMP2_y', 'CTDPRS_y', 'CTDCOND1_y', 'CTDCOND2_y', 'CTDSAL_y',
       'CTDOXY1_y', 'CTDOXYVOLTS_y', 'FREE1_y', 'FREE2_y', 'FREE3_y',
       'FREE4_y', 'FLUOR_y', 'CTDBACKSCATTER_y', 'CTDXMISS_y', 'ALT_y',
       'REF_PAR_y', 'GPSLAT_y', 'GPSLON_y', 'new_fix_y', 'pressure_temp_int_y',
       'pump_on_y', 'btl_fire_y', 'scan_datetime_y', 'SSSCC_y',
       'master_index_y','OS_y','dv_dt_y','CTDOXY',r_param]
    
    merge_df = merge_df[keep_columns]
    
    col_names = ['CTDTMP1','CTDTMP2', 'CTDPRS', 'CTDCOND1', 'CTDCOND2','CTDSAL',
       'CTDOXY1', 'CTDOXYVOLTS', 'FREE1', 'FREE2', 'FREE3',
       'FREE4', 'FLUOR', 'CTDBACKSCATTER', 'CTDXMISS', 'ALT',
       'REF_PAR', 'GPSLAT', 'GPSLON', 'new_fix', 'pressure_temp_int',
       'pump_on', 'btl_fire', 'scan_datetime', 'SSSCC',
       'master_index','OS','dv_dt','CTDOXY', r_param]
    
    merge_df.columns = [col_names]

    
    return merge_df

def get_sbe_coef(hexfile,xmlfile):


    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1, 
                        'ranking': 5, 'column': 'CTDOXYVOLTS', 
                        'sensor_info': rawConfig['Sensors'][i]}

     
    Soc = oxy_meta['sensor_info']['Soc'] # Oxygen slope
    Voff = oxy_meta['sensor_info']['offset'] # Sensor output offset voltage
    Tau20 = oxy_meta['sensor_info']['Tau20'] #Sensor time constant at 20 deg C and 1 Atm
    Tcorr = oxy_meta['sensor_info']['Tcor'] # Temperature correction
    E  = oxy_meta['sensor_info']['E'] #Compensation Coef for pressure effect on membrane permeability (Atkinson et al)
    A = oxy_meta['sensor_info']['A'] # Compensation Coef for temp effect on mem perm
    B = oxy_meta['sensor_info']['B'] # Compensation Coef for temp effect on mem perm
    C = oxy_meta['sensor_info']['C'] # Compensation Coef for temp effect on mem perm
    D = [oxy_meta['sensor_info']['D0'],oxy_meta['sensor_info']['D1'],oxy_meta['sensor_info']['D2']] # Press effect on time constant Usually not fitted for.
    
    #coef = {'Soc':Soc,'Voff':Voff,'Tau20':Tau20,'Tcorr':Tcorr,'E':E,'A':A,'B':B,'C':C,'D':D}
    
    # Make into an array in order to do keast squares fitting easier
#        coef0s:
#    coef[0] = Soc
#    coef[1] = Voffset
#    coef[2] = Tau20
#    coef[3] = Tcorr
#    coef[4] = E
#
#    cc[0] = D1
#    cc[1] = D2
    coef = [Soc,Voff,Tau20,Tcorr,E]#,A,B,C,D]
    
    return coef

def calculate_weights(pressure):
    
    eps=1e-5
    wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
    wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]#[20,20,25,25,50,50,100,100,200,200,500,500]
    wgt = scipy.interpolate.interp1d(wrow1,wrow2)
    
    weights = wgt(pressure)
    
    return weights

def interpolate_sigma(ctd_sigma,btl_sigma):
    
#    all_sigma0 = time_data_clean['sigma0_time'].append(btl_data_clean['sigma0_btl'])
#    fsg = pd.Series(np.min(all_sigma0)-1e-4)
#    lsg = pd.Series(np.max(all_sigma0)+1e-4)
#    new_sigma0 = fsg.append(time_data_clean['sigma0_time'])
#    new_sigma0 = new_sigma0.append(lsg)
    
    all_sigma = ctd_sigma.append(btl_sigma)
    fsg = pd.Series(np.min(all_sigma)-1e-4)
    lsg = pd.Series(np.max(all_sigma)+1e-4)
    new_sigma = fsg.append(ctd_sigma)
    new_sigma = new_sigma.append(lsg)
    new_sigma = new_sigma.reset_index()
    new_sigma = new_sigma[0]
    
    
    #x = np.arange(np.size(ctd_sigma))
#    x_inter = np.arange(np.size(new_sigma0))
#    inter_sigma2 = scipy.interpolate.interp1d(x_inter,new_sigma0)
#    new_x = np.linspace(0,np.max(x_inter),np.size(ctd_sigma))
#    new_sigma0=inter_sigma2(new_x)
#    new_sigma0 = pd.Series(new_sigma0)    
    
    return new_sigma

def interpolate_pressure(ctd_pres,btl_pres):
    
    all_press = ctd_pres.append(btl_pres)
    fsg = pd.Series(np.min(all_press) - 1)
    lsg = pd.Series(np.max(all_press) + 1)
    new_pres = fsg.append(ctd_pres)
    new_pres = new_pres.append(lsg)
    new_pres = new_pres.reset_index()
    new_pres = new_pres[0]
    
    return new_pres

def interpolate_param(param):
    """
    First extends 
    Generates function using x_inter
    calculates values over a new range "new_x"
    x_inter is the x-range of the extended parameter
    new_x is the x-range of the "interpolated-down" parameter 
    
    """
    ex_param = pd.Series(param.iloc[0])
    ex_param = ex_param.append(param)#,index=[])
    ex_param = ex_param.append(pd.Series(param.iloc[-1]))
    ex_param = ex_param.reset_index()
    ex_param = ex_param[0]

#    x_inter = np.arange(len(ex_param))
#    
#    inter_param = scipy.interpolate.interp1d(x_inter,ex_param)
#    
#    new_x = np.arange(len(param))
#    dparam = inter_param(new_x)
    
    return ex_param


def least_squares_resid(coef0,oxyvolts,pressure,temp,dvdt,os,ref_oxy,switch,cc=[1.92634e-4,-4.64803e-2]):
    
    coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef0,
                                     args=(oxyvolts,pressure,temp,dvdt,os,ref_oxy,switch))
    return coef

def oxygen_cal_ml(coef,oxyvolts,pressure,temp,dvdt,os,ref_oxy,switch,cc=[1.92634e-4,-4.64803e-2]):#temp,dvdt,os,

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

#    cc = [1.92634e-4,-4.64803e-2]


    #MODIFIED CODE
    
    
#    time_data['NOAA_oxy_mlL'] = coef0[0] * (time_data[oxyvo_col] 
#            + coef0[1] + coef0[2] * np.exp(cc[0] * time_data[p_col] \
#            + cc[0] * time_data[t_col]) * time_data[dvdt_col]) \
#            * time_data[os_col] * np.exp(coef0[3] * time_data[t_col]) \
#            * np.exp((coef0[4] * time_data[p_col]) \
#            / (time_data[t_col] + 273.15))
    
#    oxy_mlL = coef[0] * (oxyvolts + coef[1] + coef[2] \
#            * np.exp(cc[0] * pressure + cc[0] \
#            * temp) \
#            * dvdt) * os \
#            * np.exp(coef[3] * temp) \
#            * np.exp((coef[4] * pressure) \
#            / (temp+273.15))
#    
#    ctd_oxy_mlL = apply_oxy_coef(time_data,coef,cc=[1.92634e-4,-4.64803e-2],oxyvo_col='CTDOXYVOLTS',p_col='CTDPRS',t_col='CTDTMP1',sal_col='CTDSAL')
    ctd_oxy_mlL = oxy_from_ctd_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc)
#    ctd_oxy_mlL = SB_oxy_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc)
    #Weight Determination
    if switch == 1:
        
#        eps = 1e-5
#        
#        wrow1 = [0, 100, 100 + eps, 300, 300 + eps, 500, 500 + eps, 1200, 1200 
#                 + eps, 2000, 2000 + eps, 7000]
#        
#        wrow2 = [20, 20, 25, 25, 50, 50, 100, 100, 200, 200, 500, 500]
#        
#        wgt = scipy.interpolate.interp1d(wrow1,wrow2)
#
#        #Modified CODE
#        weights = wgt(time_data[p_col])
        
        weights = calculate_weights(pressure)

        #resid = ((weights * (ref_oxy - ctd_oxy_mlL))**2) / (np.sum(weights)**2) #Original way (np.sum(weights)**2)
        resid = ((weights * (ref_oxy - ctd_oxy_mlL))**2) #/ #(np.sum(weights)**2) #Original way (np.sum(weights)**2)
    elif switch == 2:
        #L2 Norm
        resid = (ref_oxy - ctd_oxy_mlL)**2
    
    elif switch == 3:
        #ODF residuals      
        
        resid = np.sqrt(((ref_oxy - ctd_oxy_mlL)**2) / (np.std(ctd_oxy_mlL)**2))
        
    elif switch == 4:
        # Weighted ODF residuals
        
        weights = calculate_weights(pressure)
        resid = np.sqrt(weights * ((ref_oxy - ctd_oxy_mlL)**2) / (np.sum(weights)**2))#(np.std(ctd_oxy_mlL)**2))
        
    elif switch == 5:
        
        weights = calculate_weights(pressure)

        resid = ((weights * (ref_oxy - ctd_oxy_mlL))**2) / (np.sum(weights**2))
        
        
    return resid    

def oxy_ml_to_umolkg(oxy_mlL,sigma0):
    
    oxy_uMolkg = oxy_mlL * 44660 / (sigma0 + 1000)
    
    return oxy_uMolkg


def oxy_fit(time_data, btl_data, ssscc, hexfile, xmlfile, method = 1, 
            sal_col='CTDSAL', t_col='CTDTMP1', p_col='CTDPRS', 
            lon_col='GPSLON', lat_col='GPSLAT', dov_col='CTDOXYVOLTS', 
            cond_col='CTDCOND1', sal_btl_col='CTDSAL', t_btl_col='CTDTMP1',
            p_btl_col='CTDPRS', lat_btl_col='GPSLAT', lon_btl_col='GPSLON',
            dov_btl_col='CTDOXYVOLTS', oxy_btl_col='CTDOXY1',
            cond_btl_col='CTDCOND1'):
    """
    Fits CTD oxygen data to bottle oxygen data and adds to output dataframe.
    
    Default arguments:
        
        time_data -- CTD downcast data (Pandas DataFrame)
        
        btl_data -- upcast bottle data (Pandas DataFrame)
        
        hexfile -- path to the station/cast hex file (string) 
        
        xmlfile -- path to the station/cast xml file (string)
        
    Keyword arguments:
        
        method -- least squares residual calculation method (int):
            1. Weighted L2 Norm. (Default)
            2. Vanilla L2 Norm.
            3. Normalized L2 Norm.
            
        sal_col -- DF column name for CTD downcast salinity in PSU
            (default 'CTDSAL')
        
        t_col -- DF column name for CTD downcast temperature in celcius 
            (default 'CTDTMP1')
            
        p_col -- DF column name for CTD downcast pressure in decibar
            (default 'CTDPRS')
            
        lon_col -- DF column name for CTD downcast longitude (default 'GPSLON')
        
        lat_col -- DF column name for CTD downcast latitude (default 'GPSLAT')
        
        dov_col -- DF column name for CTD downcast dissolved oxygen voltage
            (default 'CTDOXYVOLTS')
        
        cond_col -- DF column name for CTD downcast conductivity in mS/cm
            (default 'CTDCOND1')
        
        sal_btl_col -- DF column name for bottle salinity in PSU
            (default 'CTDSAL')
        
        t_btl_col -- DF column name for bottle temperature in celcius 
            (default 'CTDTMP1')
            
        p_btl_col -- DF column name for bottle pressure in decibar
            (default 'CTDPRS')
            
        lon_btl_col -- DF column name for bottle longitude (default 'GPSLON')
        
        lat_btl_col -- DF column name for bottle latitude (default 'GPSLAT')
        
        dov_btl_col -- DF column name for bottle dissolved oxygen voltage
            (default 'CTDOXYVOLTS')
        
        cond_btl_col -- DF column name for bottle conductivity in mS/cm
            (default 'CTDCOND1')
        
            
    
    """
    
    
    ##NORMALLY INT#####
    stn_nbr = int(ssscc[0:3])
    cst_nbr = int(ssscc[-2:])
    
    #Calculate Practical salinity from conductivity and put into Dataframe
    time_data[sal_col] = gsw.SP_from_C(time_data[cond_col],time_data[t_col],
                                       time_data[p_col])
    
    btl_data[sal_btl_col] = gsw.SP_from_C(btl_data[cond_col],
                                          btl_data[t_btl_col],
                                          btl_data[p_btl_col])
    
    
#    Calculate Conserv Temp, Absolute Sal, and Sigma0 for CTD data and clean up
#    the dataframe
    
    time_data['CT_time'] = gsw.CT_from_t(time_data[sal_col],time_data[t_col],
                                         time_data[p_col])
    
    time_data['SA_time'] = gsw.SA_from_SP(time_data[sal_col],time_data[p_col],
                                          time_data[lon_col],time_data[lat_col])
    
    time_data['sigma0_time'] = gsw.sigma0(time_data['SA_time'],
                                         time_data['CT_time'])
    
    time_data_clean = process_ctd.binning_df(time_data)
    time_data_clean = time_data_clean.dropna(how='all')
    time_data_clean.index = pd.RangeIndex(len(time_data_clean.index))#reindex
    
#   Rename Dataframe columns to avoid confusion later on
    time_data_clean['CTDOXYVOLTS_time'] = time_data_clean[dov_col] 
    time_data_clean['TMP_KELVIN_CTD'] = time_data_clean[t_col] + 273.15
    time_data_clean['PRESSURE_CTD'] = time_data_clean[p_col]
    time_data_clean['TEMPERATURE_CTD'] = time_data_clean[t_col]
    time_data_clean['SALINITY_CTD'] = time_data_clean[sal_col]
    
#   Reorder by sigma0    
    time_data_clean = time_data_clean.sort_values(by='sigma0_time')
    
#    Calculate Conserv Temp, Absolute Sal, and Sigma0 for BTL data and clean up
#    the dataframe    
    btl_data_clean = btl_data.dropna(how='all')
    btl_data_clean['BTLNBR'] = btl_data['btl_fire_num'] 
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index)) #reindex
    
    btl_data_clean['CT_btl'] = gsw.CT_from_t(btl_data_clean[sal_btl_col],
                                             btl_data_clean[t_btl_col],
                                             btl_data_clean[p_btl_col])
    
    btl_data_clean['SA_btl'] = gsw.SA_from_SP(btl_data_clean[sal_btl_col],
                                              btl_data_clean[p_btl_col],
                                              btl_data_clean[lon_btl_col],
                                              btl_data_clean[lat_btl_col])
    
    btl_data_clean['sigma0_btl'] = gsw.sigma0(btl_data_clean['SA_btl'],
                                              btl_data_clean['CT_btl'])
    
#   Sort bottles by sigma0
    btl_data_clean = btl_data_clean.sort_values(by='sigma0_btl') 
    
#   Calculate dV/dT for bottle data and convolve
    doxyv_btl = np.diff(btl_data_clean['CTDOXYVOLTS'])
    dt_btl = np.diff(btl_data_clean['scan_datetime'])
    dv_dt_btl = doxyv_btl/dt_btl
    dv_dt_conv_btl = np.convolve(dv_dt_btl,[0.5,0.5])
    btl_data_clean['dv_dt_conv_btl'] = dv_dt_conv_btl
    
#   Calculate dV/dT for CTD data, convolve, and filter
    doxyv = np.diff(time_data_clean['CTDOXYVOLTS_time'])
    dt = np.diff(time_data_clean['scan_datetime'])
    
    dv_dt = doxyv/dt
    dv_dt_conv = np.convolve(dv_dt,[0.5,0.5])
    
    a = 1
    windowsize = 5
    b = (1/windowsize)*np.ones(windowsize)
    
    #filt = scipy.signal.lfilter(b,a,dv_dt_conv)
    filtfilt = scipy.signal.filtfilt(b,a,dv_dt_conv)
    
    time_data_clean['dv_dt_time'] = filtfilt
    
#   Create Array containing all BTL and CTD sigma0 values and add increase the
#   range of the CTD sigma0 values in both directions by 1e-4 from the largest
#   and smallest simga0 values
    
    all_sigma0 = time_data_clean['sigma0_time'].append(btl_data_clean['sigma0_btl'])
    fsg = pd.Series(np.min(all_sigma0)-1e-4)
    lsg = pd.Series(np.max(all_sigma0)+1e-4)
    new_sigma0 = fsg.append(time_data_clean['sigma0_time'])
    new_sigma0 = new_sigma0.append(lsg)

#   Increase the size of the CTD press., sal.,temp., and oxyvo., to match the
#   length of new_sigma0
    
    new_pressure = pd.Series(time_data_clean['CTDPRS'].iloc[0])
    new_pressure = new_pressure.append(time_data_clean['CTDPRS'])
    new_pressure = new_pressure.append(pd.Series(time_data_clean['CTDPRS'].iloc[-1]))
    #len(new_pressure)

    new_sal = pd.Series(time_data_clean['CTDSAL'].iloc[0])
    new_sal = new_sal.append(time_data_clean['CTDSAL'])
    new_sal = new_sal.append(pd.Series(time_data_clean['CTDSAL'].iloc[-1]))
    #np.size(new_sal)

    new_temp = pd.Series(time_data_clean[t_col].iloc[0])
    new_temp = new_temp.append(time_data_clean[t_col])
    new_temp = new_temp.append(pd.Series(time_data_clean[t_col].iloc[-1]))
    #np.size(new_temp)

    new_oxyvo = pd.Series(time_data_clean['CTDOXYVOLTS_time'].iloc[0])
    new_oxyvo = new_oxyvo.append(time_data_clean['CTDOXYVOLTS_time'])
    new_oxyvo = new_oxyvo.append(pd.Series(time_data_clean['CTDOXYVOLTS_time'].iloc[-1]))
    #np.size(new_oxyvo)
    
#   Interpolate sigma0 for CTD Data 
    x = np.arange(np.size(time_data_clean['sigma0_time']))
    x_inter = np.arange(np.size(new_sigma0))
    inter_sigma2 = scipy.interpolate.interp1d(x_inter,new_sigma0)
    new_x = np.linspace(0,np.max(x_inter),np.size(time_data_clean['sigma0_time']))
    new_sigma0=inter_sigma2(new_x)
    
#   Interpolate press, sal, temp, and oxyvo values

    inter_pres2 = scipy.interpolate.interp1d(x_inter,new_pressure)
    inter_sal2 = scipy.interpolate.interp1d(x_inter,new_sal)
    inter_temp2 = scipy.interpolate.interp1d(x_inter,new_temp)
    inter_oxyvo2 = scipy.interpolate.interp1d(x_inter,new_oxyvo)
       
    dpres = inter_pres2(new_x)
    dsal = inter_sal2(new_x)
    dtemp = inter_temp2(new_x)
    doxyvo = inter_oxyvo2(new_x)
    
#   Add interpolated values to dataframe

    time_data_clean['sigma0_time'] = new_sigma0
    time_data_clean['PRESSURE_CTD'] = dpres
    time_data_clean['SALINITY_CTD'] = dsal
    time_data_clean['TEMPERATURE_CTD'] = dtemp
    time_data_clean['CTDOXYVOLTS_time'] = doxyvo

#   Interpolate dV/dT down to length of the bottle size
    btl_len_x = np.arange(np.size(btl_data_clean['dv_dt_conv_btl']))
    dv_dt_inter_x = np.arange(np.size(time_data_clean['dv_dt_time']))
    dv_dt_inter = np.interp(btl_len_x,dv_dt_inter_x,time_data_clean['dv_dt_time'])
    
#   Recalculate Simga0
    btl_data_clean = btl_data.dropna(how='all')
    btl_data_clean['BTLNBR'] = btl_data['btl_fire_num']
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index))
   
    btl_data_clean['CT_btl'] = gsw.CT_from_t(btl_data_clean[sal_btl_col],
                                             btl_data_clean[t_btl_col],
                                             btl_data_clean[p_btl_col])
    
    btl_data_clean['SA_btl'] = gsw.SA_from_SP(btl_data_clean[sal_btl_col],
                                              btl_data_clean[p_btl_col],
                                              btl_data_clean[lon_btl_col],
                                              btl_data_clean[lat_btl_col])
    
    btl_data_clean['sigma0_btl'] = gsw.sigma0(btl_data_clean['SA_btl'],
                                              btl_data_clean['CT_btl'])

#   Sort bottles for matching
    btl_data_clean = btl_data_clean.sort_values(by='sigma0_btl') #reindex
    
#   Merge BTL and TIME dataframes on sigma0 and reindex
    matched_df = pd.merge_asof(btl_data_clean,time_data_clean,
                               left_on='sigma0_btl',right_on='sigma0_time',
                               direction='nearest')
    
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index))

    time_data_matched = matched_df
    
#   Match CTD dV/dT direction to that of BTL data
    time_data_matched['dv_dt_time'] = dv_dt_inter*-1
    
#   Calculate Oxygen Solubility for both BTL and CTD data
    btl_data_clean['OS'] = sbe_eq.OxSol(btl_data_clean[t_btl_col],
                                        btl_data_clean['CTDSAL'])
    
    time_data_matched['OS'] = sbe_eq.OxSol(time_data_matched['TEMPERATURE_CTD'],
                                           time_data_matched['SALINITY_CTD'])

#   Determine Initial Guess from hex files    
    coef0 = None
    

    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1, 
                        'ranking': 5, 'column': 'CTDOXYVOLTS', 
                        'sensor_info': rawConfig['Sensors'][i]}

    if coef0 is None:
        coef0 = [oxy_meta['sensor_info']['Soc'], 
                 oxy_meta['sensor_info']['offset'], oxy_meta['sensor_info']['A'], 
                 oxy_meta['sensor_info']['B'], oxy_meta['sensor_info']['C'], 
                 oxy_meta['sensor_info']['E']]
        print(coef0)
        
    Tau20=oxy_meta['sensor_info']['Tau20']
    Tcorr = oxy_meta['sensor_info']['Tcor']
    
    coef0 =[oxy_meta['sensor_info']['Soc'],oxy_meta['sensor_info']['offset'],
            oxy_meta['sensor_info']['Tau20'],oxy_meta['sensor_info']['Tcor'],
            oxy_meta['sensor_info']['E']]
    
    cc=[1.92634e-4,-4.64803e-2]
    
#   Calculate oxygen in ml/L using NOAA's equation
    
    btl_data_clean['dv_dt_conv_btl'] = dv_dt_conv_btl
    
    btl_data_clean['NOAA_oxy_mlL'] = coef0[0]*(btl_data_clean['CTDOXYVOLTS'] \
        + coef0[1] + Tau20 * np.exp(cc[0] * btl_data_clean['CTDPRS'] + cc[0] 
        * btl_data_clean[t_btl_col]) \
        * btl_data_clean['dv_dt_conv_btl']) * btl_data_clean['OS'] \
        * np.exp(Tcorr * btl_data_clean[t_btl_col]) \
        * np.exp((coef0[4]*btl_data_clean['CTDPRS']) / (btl_data_clean[t_btl_col] \
        + 273.15))

    time_data_matched['NOAA_oxy_mlL'] = coef0[0] \
        * (time_data_matched['CTDOXYVOLTS_time'] \
        + coef0[1]+Tau20 * np.exp(cc[0] \
        * time_data_matched['CTDPRS_y'] + cc[0] \
        * time_data_matched['TEMPERATURE_CTD']) \
        * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
        * np.exp(Tcorr * time_data_matched['TEMPERATURE_CTD']) \
        * np.exp((coef0[4] * time_data_matched['CTDPRS_y']) \
        / (time_data_matched['TEMPERATURE_CTD'] + 273.15)) 
        
#   Calculate weights for fitting    
    eps=1e-5
    wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
    wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]
    wgt = scipy.interpolate.interp1d(wrow1,wrow2)
    
    time_data_matched['weights'] = wgt(time_data_matched['CTDPRS_y'])
    
#   Least Squares fitting to determine new coefficients   
    coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef0,
                                     args=(time_data_matched,btl_data_clean,
                                           method))

#   Apply new coeficients to recalculate Oxygen in ml/L 
    new_OXY = coef[0] * (time_data_matched['CTDOXYVOLTS_y'] + coef[1] + coef[2] \
        * np.exp(cc[0] * time_data_matched['CTDPRS_y'] + cc[0] \
        * time_data_matched['TEMPERATURE_CTD']) \
        * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
        * np.exp(coef[3] * time_data_matched['TEMPERATURE_CTD']) \
        * np.exp((coef[4] * time_data_matched['CTDPRS_y']) \
        / (time_data_matched['TEMPERATURE_CTD']+273.15))
    
#   Convert Oxygen from ml/L to uMol/kg and put CTD data into BTL dataframe
    time_data_matched['Fitted_OXY_uMOLKG_time'] = new_OXY * 44660 \
        / (time_data_matched['sigma0_time'] + 1000)
    
    btl_data_clean['OXY_uMOLKG_btl'] = btl_data_clean['NOAA_oxy_mlL'] * 44660 \
        / (btl_data_clean['sigma0_btl'] + 1000)
    
    btl_data_clean['CTDOXY'] = time_data_matched['Fitted_OXY_uMOLKG_time']
    btl_data_output = btl_data_clean
    
#   Calculate Residual and determine cutoff based upon standard deviation 
#   outlier values stored in bad_values dataframe
    btl_data_output['residual'] = np.abs(btl_data_clean['OXY_uMOLKG_btl'] 
                            - time_data_matched['Fitted_OXY_uMOLKG_time'])
 
    cutoff = np.std(btl_data_output['residual']) * 2.8
    btl_data_output['index'][btl_data_output['residual']<=cutoff]
    btl_data_clean=btl_data_output[btl_data_output['residual']<=cutoff]
    thrown_values = btl_data_output[btl_data_output['residual']>cutoff]
    bad_values = btl_data_output[btl_data_output['residual']>cutoff]
    btl_data_clean
    print(thrown_values['BTLNBR'])
    
#   Loop until there are no outliers
    
    while not thrown_values.empty:
        
#       Re-match new BTL dataframe with original time dataframe and reindex
        
        matched_df = pd.merge_asof(btl_data_clean,time_data_clean,
                                   left_on='sigma0_btl',right_on='sigma0_time',
                                   direction='nearest')
        
        btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index)) #reindex DF
        time_data_matched = matched_df
        
#       Re-Interpolate dV/dT
        btl_len_x = np.arange(np.size(btl_data_clean['dv_dt_conv_btl']))
        dv_dt_inter_x = np.arange(np.size(time_data_clean['dv_dt_time']))
        
        dv_dt_inter = np.interp(btl_len_x,dv_dt_inter_x,
                                time_data_clean['dv_dt_time'])
        
        time_data_matched['dv_dt_time'] = dv_dt_inter*-1
        
#       Recalculate Oxygen Solubility

        btl_data_clean['OS'] = sbe_eq.OxSol(btl_data_clean[t_btl_col],
                                            btl_data_clean['CTDSAL'])
        
        time_data_matched['OS'] = sbe_eq.OxSol(time_data_matched['TEMPERATURE_CTD'],
                                               time_data_matched['SALINITY_CTD'])
        
#       Recalculate Oxygen (ml/L)       
        time_data_matched['NOAA_oxy_mlL'] = coef0[0] \
                * (time_data_matched['CTDOXYVOLTS_time'] \
                + coef0[1]+Tau20 * np.exp(cc[0] \
                * time_data_matched['CTDPRS_y'] + cc[0] \
                * time_data_matched['TEMPERATURE_CTD']) \
                * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
                * np.exp(Tcorr * time_data_matched['TEMPERATURE_CTD']) \
                * np.exp((coef0[4] * time_data_matched['CTDPRS_y']) \
                / (time_data_matched['TEMPERATURE_CTD'] + 273.15))   
                
#       Recalculate coeficients 
        coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef,
                                         args=(time_data_matched,btl_data_clean,
                                               method))
        
#       Re-apply coeficients
        new_OXY = coef[0] * (time_data_matched['CTDOXYVOLTS_y'] + coef[1] + coef[2] \
                  * np.exp(cc[0] * time_data_matched['CTDPRS_y'] + cc[0] \
                  * time_data_matched['TEMPERATURE_CTD']) \
                  * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
                  * np.exp(coef[3] * time_data_matched['TEMPERATURE_CTD']) \
                  * np.exp((coef[4] * time_data_matched['CTDPRS_y']) \
                  / (time_data_matched['TEMPERATURE_CTD']+273.15))
                  
#       Convert oxygen ml/L to uMol/kg
        time_data_matched['Fitted_OXY_uMOLKG_time'] = new_OXY * 44660 \
                  / (time_data_matched['sigma0_time'] + 1000)
                  
        btl_data_clean['OXY_uMOLKG_btl'] = btl_data_clean['NOAA_oxy_mlL'] \
                  * 44660 / (btl_data_clean['sigma0_btl'] + 1000)
                  
#       Rename Columns
        btl_data_clean['OXYGEN'] = btl_data_clean['OXY_uMOLKG_btl']
        btl_data_clean['CTDOXY'] = time_data_matched['Fitted_OXY_uMOLKG_time']
        
#       Recalculate residuals and residual standard deviation
        btl_data_output = btl_data_clean
        btl_data_clean['residual'] = np.abs(btl_data_output['OXYGEN']-btl_data_output['CTDOXY'])
        std_res = np.std(btl_data_output['residual'])
        cutoff = np.std(btl_data_output['residual'])*2.8
                  
#       Use residuals to determine cutoff and collect outliers
        btl_data_clean=btl_data_clean[btl_data_output['residual']<=cutoff]
        bad_values2 = btl_data_output[btl_data_output['residual']>cutoff]
        thrown_values = btl_data_clean[btl_data_output['residual']>cutoff]
        print(thrown_values['BTLNBR'])
        print(btl_data_clean['BTLNBR'])
        bad_values = pd.concat([bad_values,bad_values2])
        print('NEW STANDARD DEVIATION IS:',std_res)
        
#      Add Station and Cast Numbers to DataFrame
    time_data_matched['STNNBR']=stn_nbr
    time_data_matched['CASTNO']=cst_nbr

    btl_data_clean['STNNBR']=stn_nbr
    btl_data_clean['CASTNO']=cst_nbr
        
#   Sanity PLOT
#    plt.plot(btl_data_clean['residual'],btl_data_clean['CTDPRS']*-1,'bx')
#    plt.xlim(xmax=10)
#    plt.show()
    #dataframe_concat = pd.concat([dataframe_concat,btl_data_clean])
        
#   Flag data    
    bad_values['CTDOXY_FLAG_W'] = 4
    btl_data_clean['CTDOXY_FLAG_W'] = 2
    btl_data_clean = pd.concat([btl_data_clean,bad_values])
    btl_data_clean=btl_data_clean.sort_values(by='BTLNBR')
    
#   Compose writeable dataframe for exportation 
    btl_data_clean = btl_data_clean.reset_index(drop=True)
    btl_data_write = pd.DataFrame()
    btl_data_write['STNNBR'] = btl_data_clean['STNNBR'].values
    btl_data_write['CASTNO'] = btl_data_clean['CASTNO'].values
    btl_data_write['SAMPNO'] = btl_data_clean['BTLNBR'].values
    btl_data_write['CTDOXY'] = btl_data_clean['CTDOXY'].values
    btl_data_write['CTDOXY_FLAG_W'] = btl_data_clean['CTDOXY_FLAG_W'].values
    
#    btl_data_write.to_csv(filestring,index=False)
#    btl_write_concat = pd.concat([btl_write_concat,btl_data_write])
    
    return coef #btl_data_write, btl_data_clean, coef


    
    
    
    #####           GRAVEYARD             ##########
        
    #NOAA Calucaltion (unmodified):
        #ORIGINAL CODE
        #time_data['NOAA_oxy_mlL'] =coef0[0]*(time_data['CTDOXYVOLTS_y']+\
        #coef0[1]+coef0[2]*np.exp(cc[0]*time_data['CTDPRS_y']+cc[0]*time_data['TEMPERATURE_CTD'])\
        #*time_data['dv_dt_time'])*time_data['OS']*np.exp(coef0[3]*time_data['TEMPERATURE_CTD'])\
        #*np.exp((coef0[4]*time_data['CTDPRS_y'])/(time_data['TEMPERATURE_CTD']+273.15))
    
    #SWITCH 1 code:
    
        #Original Code
        #time_data['weights'] = wgt(time_data['CTDPRS_y'])
        #resid = np.sum((time_data['weights']*(btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL']))/(np.sum(time_data['weights'])**2))
        #resid = (time_data['weights']*(btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL']))/(np.sum(time_data['weights'])**2) Working
        
    #SWITCH 3 code:
    
        #Original Code
        #resid = np.sqrt(((btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL'])**2)/(np.std(time_data['NOAA_oxy_mlL'])**2))