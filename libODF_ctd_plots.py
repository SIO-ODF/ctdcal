import matplotlib
matplotlib.use('agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#####
# The following section is a little brittle due to hardcoded names, but we'll fix
# that later. Code copy pasted from jupyter notebook.
#####

def all_plots(df):
    '''Create and output all plots'''
    btl_t1_residuals_pressure_plot(df)
    btl_t2_residuals_pressure_plot(df)
    t1_t2_residuals_pressure_plot(df)
    btl_t1_residuals_station_plot(df)
    btl_t2_residuals_station_plot(df)
    t1_t2_residuals_station_plot(df)
    btl_t1_residuals_station_deep_plot(df)
    btl_t2_residuals_station_deep_plot(df)
    t1_t2_residuals_station_deep_plot(df)
    btl_c1_residuals_pressure_plot(df)
    btl_c2_residuals_pressure_plot(df)
    c1_c2_residuals_pressure_plot(df)
    btl_c1_residuals_station_plot(df)
    btl_c2_residuals_station_plot(df)
    c1_c2_residuals_station_plot(df)
    btl_c1_residuals_station_deep_plot(df)
    btl_c2_residuals_station_deep_plot(df)
    c1_c2_residuals_station_deep_plot(df)
    c_t_coherence_plot(df)
    btl_c1_residuals_compare_plot(df)
    btl_c2_residuals_compare_plot(df)
    c1_c2_residuals_compare_plot(df)
    btl_c1_residuals_station_uncorrected_plot(df)
    btl_c2_residuals_station_uncorrected_plot(df)
    c1_c2_residuals_station_uncorrected_plot(df)
    btl_sal_pressure_plot(df)
    btl_sal_station_plot(df)
    btl_sal_station_deep_plot(df)
    btl_oxy_residuals_pressure_plot(df)
    btl_oxy_residuals_station_plot(df)
    btl_oxy_residuals_station_deep_plot(df)
    return None

     #################################################################
     ##### Here lies the temperature plots, long may they rest.  #####
     #################################################################

def btl_t1_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_T1'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('REFTMP-CTDTMP1 vs CTDPRS')
    ax.set_xlabel('T1 Residual (T90 C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/reftmp_t1_p.ps', format='ps')
    fig.savefig('./images/reftmp_t1_p.png', format='png')
    plt.close()
    return None

def btl_t2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_T2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('REFTMP-CTDTMP2 vs CTDPRS')
    ax.set_xlabel('T2 Residual (T90 C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/reftmp_t2_p.ps', format='ps')
    fig.savefig('./images/reftmp_t2_p.png', format='png')
    plt.close()
    return None

def t1_t2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['T1_T2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('CTDTMP1-CTDTMP2 vs CTDPRS')
    ax.set_xlabel('T1-T2 Residual (T90 C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/t1_t2_p.ps', format='ps')
    ig.savefig('./images/t1_t2_p.png', format='png')
    plt.close()
    return None

def btl_t1_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_T1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP1 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1 Residual (T90 C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/reftmp_t1_stn.ps', format='ps')
    fig.savefig('./images/reftmp_t1_stn.png', format='png')
    plt.close()
    return None

def btl_t2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_T2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T2 Residual (T90 C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/reftmp_t2_stn.ps', format='ps')
    fig.savefig('./images/reftmp_t2_stn.png', format='png')
    plt.close()
    return None

def t1_t2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_T1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDTMP1-CTDTMP2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1-T2 Residual (T90 C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/t1_t2_stn.ps', format='ps')
    fig.savefig('./images/t1_t2_stn.png', format='png')
    plt.close()
    return None

def btl_t1_residuals_station_deep_plot(df):
    '''Decide how to deal with deepness - should the dataframe be already cut to
    2000 dbar or not? If variable, it should be a kwargs. Update comments in others
    df in must be cut to all samples > 2000 dbar
    Ex: df = df[df['CTDPRS'] > 2000]
    '''
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['BTL_T1'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP1 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1 Residual (T90 C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/reftmp_t1_stn_deep.ps', format='ps')
    fig.savefig('./images/reftmp_t1_stn_deep.png', format='png')
    plt.close()
    return None

def btl_t2_residuals_station_deep_plot(df):
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['BTL_T2'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T2 Residual (T90 C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/reftmp_t2_stn_deep.ps', format='ps')
    fig.savefig('./images/reftmp_t2_stn_deep.png', format='png')
    plt.close()
    return None

def t1_t2_residuals_station_deep_plot(df):
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['T1_T2'], marker='+', c=df_deep['CTDPRS'], cmap='viridis_r')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDTMP1-CTDTMP2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1-T2 Residual (T90 C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/t1_t2_stn_deep.ps', format='ps')
    fig.savefig('./images/t1_t2_stn_deep.png', format='png')
    plt.close()
    return None

     #################################################################
     ##### Here lies the conductivity plots, long may they rest. #####
     #################################################################

def btl_c1_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_C1'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('BTLCOND-CTDCOND1 vs CTDPRS')
    ax.set_xlabel('C1 Residual (mS/cm)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/btlcond_c1_p.ps', format='ps')
    fig.savefig('./images/btlcond_c1_p.png', format='png')
    plt.close()
    return None

def btl_c2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_C2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('BTLCOND-CTDCOND2 vs CTDPRS')
    ax.set_xlabel('C2 Residual (mS/cm)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/btlcond_c2_p.ps', format='ps')
    fig.savefig('./images/btlcond_c2_p.png', format='png')
    plt.close()
    return None

def c1_c2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['C1_C2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('CTDCOND1-CTDCOND2 vs CTDPRS')
    ax.set_xlabel('C1-C2 Residual (mS/cm)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/c1_c2_p.ps', format='ps')
    fig.savefig('./images/c1_c2_p.png', format='png')
    plt.close()
    return None

def btl_c1_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND-CTDCOND1 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c1_stn.ps', format='ps')
    fig.savefig('./images/btlcond_c1_stn.png', format='png')
    plt.close()
    return None

def btl_c2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND-CTDCOND2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c2_stn.ps', format='ps')
    fig.savefig('./images/btlcond_c2_stn.png', format='png')
    plt.close()
    return None

def c1_c2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['C1_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDCOND1-CTDCOND2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1-C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/c1_c2_stn.ps', format='ps')
    fig.savefig('./images/c1_c2_stn.png', format='png')
    plt.close()
    return None

def btl_c1_residuals_station_deep_plot(df):
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['BTL_C1'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND-CTDCOND1 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c1_stn_deep.ps', format='ps')
    fig.savefig('./images/btlcond_c1_stn_deep.png', format='png')
    plt.close()
    return None

def btl_c2_residuals_station_deep_plot(df):
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['BTL_C2'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND-CTDCOND2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c2_stn_deep.ps', format='ps')
    fig.savefig('./images/btlcond_c2_stn_deep.png', format='png')
    plt.close()
    return None

def c1_c2_residuals_station_deep_plot(df):
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['C1_C2'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDCOND1-CTDCOND2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1-C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/c1_c2_stn_deep.ps', format='ps')
    fig.savefig('./images/c1_c2_stn_deep.png', format='png')
    plt.close()
    return None

def c_t_coherence_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['T1_T2'],df['C1_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('T1-T2 vs C1-C2')
    ax.set_xlabel('T1-T2 Residual (T90 C)')
    ax.set_ylabel('C1-C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/c_t_coherence_p.ps', format='ps')
    fig.savefig('./images/c_t_coherence_p.png', format='png')
    plt.close()
    return None

def btl_c1_residuals_compare_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTLCOND'], df['BTL_C1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND vs BTLCOND-CTDCOND1')
    ax.set_xlabel('Reference Conductivity (mS/cm)')
    ax.set_ylabel('C1 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c1_compare.ps', format='ps')
    fig.savefig('./images/btlcond_c1_compare.png', format='png')
    plt.close()
    return None

def btl_c2_residuals_compare_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTLCOND'], df['BTL_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND vs BTLCOND-CTDCOND2')
    ax.set_xlabel('Reference Conductivity (mS/cm)')
    ax.set_ylabel('C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c2_compare.ps', format='ps')
    fig.savefig('./images/btlcond_c2_compare.png', format='png')
    plt.close()
    return None

def c1_c2_residuals_compare_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTLCOND'], df['C1_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND vs CTDCOND1-CTDCOND2')
    ax.set_xlabel('Reference Conductivity (mS/cm)')
    ax.set_ylabel('C1-C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/c1_c2_compare.ps', format='ps')
    fig.savefig('./images/c1_c2_compare.png', format='png')
    plt.close()
    return None

def btl_c1_residuals_station_uncorrected_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND-CTDCOND1 (Uncorrected) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c1_stn_uncorrected.ps', format='ps')
    fig.savefig('./images/btlcond_c1_stn_uncorrected.png', format='png')
    plt.close()
    return None

def btl_c2_residuals_station_uncorrected_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('BTLCOND-CTDCOND2 (Uncorrected) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlcond_c2_stn_uncorrected.ps', format='ps')
    fig.savefig('./images/btlcond_c2_stn_uncorrected.png', format='png')
    plt.close()
    return None

def c1_c2_residuals_station_uncorrected_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['C1_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDCOND1-CTDCOND2 (Uncorrected) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1-C2 Residual (mS/cm)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/c1_c2_stn_uncorrected.ps', format='ps')
    fig.savefig('./images/c1_c2_stn_uncorrected.png', format='png')
    plt.close()
    return None

def btl_sal_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_SAL'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('SALNTY-CTDSAL vs CTDPRS')
    ax.set_xlabel('CTDSAL Residual (mPSU)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/btlsal_sal_p.ps', format='ps')
    fig.savefig('./images/btlsal_sal_p.png', format='png')
    plt.close()
    return None

def btl_sal_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_SAL'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('SALNTY-CTDSAL vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('CTDSAL Residual (mPSU)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlsal_sal_stn.ps', format='ps')
    fig.savefig('./images/btlsal_sal_stn.png', format='png')
    plt.close()
    return None

def btl_sal_station_deep_plot(df):
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['BTL_SAL'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('SALNTY-CTDSAL (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('CTDSAL Residual (mPSU)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btlsal_sal_stn_deep.ps', format='ps')
    fig.savefig('./images/btlsal_sal_stn_deep.png', format='png')
    plt.close()
    return None

     #################################################################
     ######## Here lies the oxygen plots, long may they rest. ########
     #################################################################

def btl_oxy_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_O'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-10,10)
    ax.set_title('OXYGEN-CTDOXY vs CTDPRS')
    ax.set_xlabel('CTDOXY Residual (umol/kg)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./images/btl_oxy_p.ps', format='ps')
    fig.savefig('./images/btl_oxy_p.png', format='png')
    plt.close()
    return None

def btl_oxy_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'],df['BTL_O'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-10,10)
    ax.set_title('OXYGEN-CTDOXY vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('CTDOXY Residual (umol/kg)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btl_oxy_stn.ps', format='ps')
    fig.savefig('./images/btl_oxy_stn.png', format='png')
    plt.close()
    return None

def btl_oxy_residuals_station_deep_plot(df):
    df_deep = df[df['CTDPRS'] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'],df_deep['BTL_O'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-10,10)
    ax.set_title('OXYGEN-CTDOXY (> 2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('CTDOXY Residual (umol/kg)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')

    fig.savefig('./images/btl_oxy_stn_deep.ps', format='ps')
    fig.savefig('./images/btl_oxy_stn_deep.png', format='png')
    plt.close()
    return None
