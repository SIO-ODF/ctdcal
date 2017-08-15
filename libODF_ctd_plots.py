import matplotlib
matplotlib.use('nbagg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#####
# The following section is a little brittle due to hardcoded names, but we'll fix
# that later. Code copy pasted from jupyter notebook.
#####

     #################################################################
     ##### Here lies the temperature plots, long may they rest.  #####
     #################################################################

def btl_t1_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_T1'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('REFTMP-CTDTMP1 vs CTDPRS')
    ax.set_xlabel('T1 Residual (T90 milliDegrees C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')

    fig.savefig('./reftmp_t1_p.ps', format='ps')
    plt.close()
    return None

def btl_t2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_T2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('REFTMP-CTDTMP2 vs CTDPRS')
    ax.set_xlabel('T2 Residual (T90 milliDegrees C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')
    fig.savefig('./reftmp_t2_p.ps', format='ps')
    plt.close()
    return None

def t1_t2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['T1_T2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('CTDTMP1-CTDTMP2 vs CTDPRS')
    ax.set_xlabel('T1-T2 Residual (T90 milliDegrees C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')
    fig.savefig('./t1_t2_p.ps', format='ps')
    plt.close()
    return None

def btl_t1_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_T1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP1 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./reftmp_t1_stn.ps', format='ps')
    plt.close()
    return None

def btl_t2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_T2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./reftmp_t2_stn.ps', format='ps')
    plt.close()
    return None

def t1_t2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_T1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDTMP1-CTDTMP2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1-T2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./t1_t2_stn.ps', format='ps')
    plt.close()
    return None

def btl_t1_residuals_station_deep_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['BTL_T1'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP1 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./reftmp_t1_stn_deep.ps', format='ps')
    plt.close()
    return None

def btl_t2_residuals_station_deep_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['BTL_T2'], marker='+', c=df_deep['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('REFTMP-CTDTMP2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./reftmp_t2_stn_deep.ps', format='ps')
    plt.close()
    return None

def t1_t2_residuals_station_deep_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df_deep['STNNBR'], df_deep['T1_T2'], marker='+', c=df_deep['CTDPRS'], cmap='viridis_r')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDTMP1-CTDTMP2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('T1-T2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./t1_t2_stn_deep.ps', format='ps')
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
    ax.set_title('SALNTY-CTDCOND1 vs CTDPRS')
    ax.set_xlabel('C1 Residual (T90 milliDegrees C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')
    fig.savefig('./btlcond_c1_p.ps', format='ps')
    plt.close()
    return None

def btl_c2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['BTL_C2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('SALNTY-CTDCOND2 vs CTDPRS')
    ax.set_xlabel('C2 Residual (T90 milliDegrees C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')
    fig.savefig('./btlcond_c2_p.ps', format='ps')
    plt.close()
    return None

def c1_c2_residuals_pressure_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['C1_C2'],-df['CTDPRS'], marker='+', c=df['STNNBR'], cmap='rainbow')
    ax.set_xlim(-0.02,0.02)
    ax.set_title('CTDCOND1-CTDCOND2 vs CTDPRS')
    ax.set_xlabel('C1-C2 Residual (T90 milliDegrees C)')
    ax.set_ylabel('Pressure (dbar)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Station Number')
    fig.savefig('./c1_c2_p.ps', format='ps')

def btl_c1_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('SALNTY-CTDCOND1 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./btlcond_c1_stn.ps', format='ps')
    plt.close()
    return None

def btl_c2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('SALNTY-CTDCOND2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./btlcond_c2_stn.ps', format='ps')
    plt.close()
    return None

def c1_c2_residuals_station_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['C1_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDCOND1-CTDCOND2 vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1-C2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./c1_c2_stn.ps', format='ps')
    plt.close()
    return None

def btl_c1_residuals_station_deep_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C1'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('SALNTY-CTDCOND1 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./btlcond_c1_stn_deep.ps', format='ps')
    plt.close()
    return None

def btl_c2_residuals_station_deep_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['BTL_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('SALNTY-CTDCOND2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./btlcond_c2_stn_deep.ps', format='ps')
    plt.close()
    return None

def c1_c2_residuals_station_deep_plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    cm = ax.scatter(df['STNNBR'], df['C1_C2'], marker='+', c=df['CTDPRS'], cmap='rainbow')
    ax.set_ylim(-0.01,0.01)
    ax.set_title('CTDCOND1-CTDCOND2 (>2000 db) vs STNNBR')
    ax.set_xlabel('Station Number')
    ax.set_ylabel('C1-C2 Residual (T90 milliDegrees C)')
    cbar = fig.colorbar(cm)
    cbar.set_label('Pressure (dbar)')
    fig.savefig('./c1_c2_stn_deep.ps', format='ps')
    plt.close()
    return None
