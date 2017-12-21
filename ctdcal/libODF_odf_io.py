import pandas as pd
import numpy as np
import datetime as dt

import sys

def salts_time_indexer(time_string):
    '''Take in string time, then change it to seconds relative to a 0 point.'''
    try:
        output = dt.datetime.strptime(time_string, f'%H:%M:%S')
    except ValueError:
        output = dt.datetime.strptime(time_string, f'%H:%M:%S*')
    output = output - dt.datetime(1900, 1, 1, 0, 0, 0)
    return output.seconds

def load_short_salts_2017(filepath):
    '''Load a salts file in truncated form and prep for pre-calibration. c. Sept 2017
    Follows normal form.

    Input:
    filepath - a string or path object to the file
    Output:
    dataframe - df following normal form
    '''

    names = [
    'STNNBR',
    'CASTNO',
    'autosal_box_num',
    'bath_temp',
    'cr_average',
    'SAMPNO',
    'standby_number',
    'end_time'
    ]
    usecols = [0,1,2,3,4,5,6,8]

    dtype = {
        'STNNBR':np.int32,
        'CASTNO':np.int32,
        'autosal_box_num':np.int32,
        'bath_temp':np.int32,
        'cr_average':np.float64,
        'SAMPNO':object,
        'standby_number':np.int32
    }
    converters = {
        'end_time': salts_time_indexer
    }

    df = pd.read_csv(filepath, usecols = usecols, names = names, skiprows = [0],
                     delim_whitespace = True, dtype = dtype, converters = converters)

    return df

def salts_create_index_time(df):
    '''Deal with runs that go through midnight, but end less than 24 hours later.'''
    df['index_time'] = df['end_time'] - df['end_time'].iloc[0]
    #make a boolean array of 0s and 1s, then multiply the array by the number to get the number you want, then add together
    df["index_time"] += (df["index_time"] < 0) * (3600*24)
    return df

def create_salt_std_time_array(df):
    '''Create a preformatted array to be passed to np.polyfit for x, or time'''
    df = df.groupby(['SAMPNO'])
    df = df.get_group('worm')
    return df['index_time']

def create_salt_std_cr_array(df):
    '''Create a preformatted array to be passed to np.polyfit for y, or cr'''
    df = df.groupby(['SAMPNO'])
    df = df.get_group('worm')
    return df['cr_average']

def autosal_drift_calibration(df):
    '''Calculate a linear polynomial fit for autosal drift between start and end of runs.'''
    p = np.polyfit(create_salt_std_time_array(df), create_salt_std_cr_array(df), 1)
    p[1] = 0
    p = np.poly1d(p)
    return p

def autosal_drift_fit(df):
    '''Take in the base io dataframe, then compute calibration, then apply fit to data'''
    df = salts_create_index_time(df)
    p = autosal_drift_calibration(df)
    df['cr_average_fitted'] = p(df['index_time']) + df['cr_average']
    df = df[df['SAMPNO'] != 'worm']
    return df

'''The following extracted from gsw library before conversion to python wrapper of C'''

def SP_salinometer(Rt, t):
    r"""Calculates Practical Salinity SP from a salinometer, primarily using
    the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical
    Salinity is only valid in the range 2 < SP < 42.  If the PSS-78 algorithm
    produces a Practical Salinity that is less than 2 then the Practical
    Salinity is recalculated with a modified form of the Hill et al. (1986)
    formula. The modification of the Hill et al. (1986) expression is to
    ensure that it is exactly consistent with PSS-78 at SP = 2.

    A laboratory salinometer has the ratio of conductivities, Rt, as an output,
    and the present function uses this conductivity ratio and the temperature t
    of the salinometer bath as the two input variables.

    Parameters
    ----------
    Rt : array
         C(SP,t_68,0)/C(SP=35,t_68,0) [unitless]
         conductivity ratio
         :math:`R = \frac{C(S, t_68, 0)}{C(35, 15(IPTS-68),0)} [unitless]

    t : array
        Temperature of the bath of the salinometer [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    SP : array
         Practical Salinity [psu (PSS-78), unitless]

    See Also
    --------
    TODO: sw.sals

    Notes
    -----
    TODO

    Examples
    --------
    TODO

    References
    -----------
    ..[1] Fofonoff, P. and R.C. Millard Jr. 1983: Algorithms for computation of
    fundamental properties of seawater.  Unesco Tech. Pap. in Mar. Sci., 44,
    53 pp.

    ..[2] Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the
    Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng., 11,
    109 - 112.

    .. [3] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See appendix E of this TEOS-10 Manual, and in
    particular, Eqns. (E.2.1) and (E.2.6).

    Modifications:
    2011-04-30. Paul Barker, Trevor McDougall and Rich Pawlowicz. Version 3.0
    """
    a = (0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081)
    b = (0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144)
    k = 0.0162

    Rt, t = np.broadcast_arrays(Rt, t)

    t68 = t * 1.00024
    ft68 = (t68 - 15) / (1 + k * (t68 - 15))

    Rt[Rt < 0] = np.ma.masked
    Rtx = np.sqrt(Rt)

    SP = (a[0] + (a[1] + (a[2] + (a[3] + (a[4] + a[5] * Rtx) * Rtx) * Rtx) *
                  Rtx) * Rtx + ft68 *
          (b[0] + (b[1] + (b[2] + (b[3] + (b[4] + b[5] * Rtx) * Rtx) * Rtx) *
                   Rtx) * Rtx))

    """The following section of the code is designed for SP < 2 based on the
    Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
    exactly equal to the PSS-78 algorithm at SP = 2."""

    I2 = SP < 2
    if I2.any():
        Hill_ratio = Hill_ratio_at_SP2(t[I2])
        x = 400 * Rt[I2]
        sqrty = 10 * Rtx[I2]
        part1 = 1 + x * (1.5 + x)
        part2 = 1 + sqrty * (1 + sqrty * (1 + sqrty))
        SP_Hill_raw = SP[I2] - a[0] / part1 - b[0] * ft68[I2] / part2
        SP[I2] = Hill_ratio * SP_Hill_raw
    # Ensure that SP is non-negative.
    SP = np.maximum(SP, 0)
    return SP

def compute_salinity(df):
    df['SALNTY'] = SP_salinometer(df['cr_average_fitted']/2, df['bath_temp'])
    return df

def formatted_salt_file(df):
    df = df[['STNNBR', 'CASTNO', 'SAMPNO', 'SALNTY']]
    return df

def main(argv):
    '''Example script on how to run functions. Not intended for use.'''
    #argv should be a filename
    df = load_short_salts_2017(f'{argv}')
    df = salts_create_index_time(df)
    p = autosal_drift_calibration(df)
    df = autosal_drift_fit(df)
    df = compute_salinity(df)
    df = formatted_salt_file(df)
    return None

if __name__ == '__main__':
    main(sys.argv[1:])
