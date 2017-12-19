import pandas as pd
import numpy as np

def simple_filter(df, **kwargs):
    '''Recreate simple filter for pressure/residual ranges.

    The filter spec is a dict of three values, ex: {'start_pressure': int, 'max_pressure': int, 'range': int}
    'start_pressure' is the lower pressure, inclusive
    'max_pressure' is the higher pressure, exclusive
    'range' is the acceptable limits for differences between bottle and ctd values
    The filter can be a dict or a list of dicts, and will apply all filters to the data.

    The parameter difference spec is a list of two values, ex: [str, str]
    'value1' and 'value2' are the names of columns that currently exist in the dataframe
    'value1_value2' column will be created in the dataframe with the difference of the two columns
    The parameter difference can be a list or a list of lists, and will compute all differences.

    Output:
    df - a dataframe that holds datapoints flagged as 3 until manual review.
    Spec:
    'STNNBR' 'CASTNO' 'SAMPNO' 'parameter' 'flag'
    '''



    pressure_filter = kwargs.get('filter', [
        {'start_pressure': 2000, 'max_pressure': 10000, 'range': 0.002},
        {'start_pressure': 1000, 'max_pressure': 2000, 'range': 0.005},
        {'start_pressure': 500, 'max_pressure': 1000, 'range': 0.010},
        {'start_pressure': 0, 'max_pressure': 500, 'range': 0.020}
        ]
        )
    param_diff = kwargs.get('diff', [
        ['REFTMP','CTDTMP1'],
        ['REFTMP','CTDTMP2'],
        ['CTDTMP1','CTDTMP2'],
        ['BTLCOND','CTDCOND1'],
        ['BTLCOND','CTDCOND2'],
        ['CTDCOND1','CTDCOND2'],
        ['SALNTY','CTDSAL']
        ]
        )

    compute_btlcond(df)
    diff_columns = []
    #compute diffs
    for parameter_pair in param_diff:
        diff_columns.append(general_diff(df, parameter_pair))
    #group by parameter in a cast instead of by pressures in a cast
    for column in diff_columns:
        for filter_values in pressure_filter:
            pressure_filter(df, filter_values, column)

    #df_flagged = df[]
    return None

def general_diff(df, param_list):
    '''Compute residuals and add to dataframe.
    Operate in place.
    '''
    df[f'{param_list[0]}_{param_list[1]}'] = df[f'{param_list[0]}'] - df[f'{param_list[1]}']
    return f'{param_list[0]}_{param_list[1]}'

def compute_btlcond(df):
    '''Compute bottle conductivities from salinity. Need to determine how to do it from '''
    try:
        df['BTLCOND'] = gsw.C_from_SP(df['SALNTY'], df['CTDTMP1'], df['CTDPRS'])
    except ValueError:
        print(f'Could not compute bottle conductivities, fixing -999/NaN values')
        df['SALNTY'] = df['SALNTY'].replace(to_replace = -999, value = np.nan)
        df['BTLCOND'] = gsw.C_from_SP(df['SALNTY'], df['CTDTMP1'], df['CTDPRS'])
        return None

def pressure_filter(df, filter_dict, column_name):
    '''Compare values and filter differences that may be too big.
    Operates on pressure, expecting pressure column to be called 'CTDPRS'.
    Only operates on one pair of parameters at a time - business logic in something above.
    Output -
    Return a dataframe formatted as ['STNNBR','CASTNO','SAMPNO',f'{column_name}',f'{column_name}_flag']
    The input dataframe will not be modified
    '''
    ### in case the flag column hasn't been created, make it now ###
    df[f'{column_name}_flag'] = 2
    #find all samples between the pressures we want
    df_view = df[(df['CTDPRS'] < filter_dict['max_pressure']) & (df['CTDPRS'] >= filter_dict['start_pressure'])]
    #find all samples whose value fall outside of the range we want
    df_view = df_view[(df_view[f'{column_name}'] > filter_dict['range']) | (df_view[f'{column_name}'] < -filter_dict['range'])]
    #code those samples 3
    df_view[f'{column_name}_flag'] = 3

    ### TODO: Additional logic may go here ###

    return df_view['STNNBR','CASTNO','SAMPNO',f'{column_name}',f'{column_name}_flag']

def pressure_filter_all(df, filter_dict, diff_columns):
    pass
