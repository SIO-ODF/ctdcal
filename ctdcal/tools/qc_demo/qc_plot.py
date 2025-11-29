#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Displays bottle data by parameter with plots for parameter vs pressure and residuals.
"""
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

from ctdcal.common import make_cast_id_list
from ctdcal.tools.qc_demo.qc_demo_app import store_persistent

# definitions
datadir = st.session_state.exchange_dir
file_list = make_cast_id_list(datadir, pattern='*ct1.csv')
cast_list = [str(fname)[:-4] for fname in file_list]

btl_file = list(datadir.glob('*hy1*'))[0]
btl_data = pd.read_csv(
        btl_file,
        header=0, skiprows=[0, 2], skipfooter=1,
        parse_dates={'datetime': ['DATE', 'TIME']},
        na_values='-999',
        engine='python',
)
btl_data['cast_id'] = (
        btl_data['STNNBR'].astype(str).str.rjust(3, '0') +
        btl_data['CASTNO'].astype(str).str.rjust(2, '0')
)

param_names = {
        'CTDPRS': 'Pressure',
        'CTDTMP': 'CTD Temperature',
        'REFTMP': 'Ref Temperature',
        'CTDSAL': 'CTD Salinity',
        'SALNTY': 'Bottle Salinity',
        'CTDOXY': 'CTD Oxygen',
        'OXYGEN': 'Bottle Oxygen',
}


# Begin drawing the page
st.set_page_config(page_title='QC Plots', page_icon=':ocean:', layout='wide')
page_left, page_right = st.columns([1, 1])

with page_left:
    # cast selection dropdown
    if 'persistent_cast' not in st.session_state:
        st.session_state.selected_cast = cast_list[0]
    else:
        st.session_state.selected_cast = st.session_state.persistent_cast

    st.selectbox(
            'Cast',
            cast_list,
            # index=cast_list.index(st.session_state.selected_cast),
            key='selected_cast',
            on_change=store_persistent,
            args=('selected_cast', 'persistent_cast'),
    )
    st.subheader('Viewing cast %s' % st.session_state.selected_cast)
    # filter the selected cast
    data = btl_data.loc[btl_data['cast_id'] == st.session_state.selected_cast]

    # set up parameter selection and draw the selection bar
    parameter_options = ['Salinity', 'Oxygen', 'Temperature']
    if 'persistent_parameter' not in st.session_state:
        st.session_state.selected_parameter = parameter_options[0]
    else:
        st.session_state.selected_parameter = st.session_state.persistent_parameter

    st.segmented_control(
            "Parameter",
            options=parameter_options,
            selection_mode="single",
            key='selected_parameter',
            on_change=store_persistent,
            args=('selected_parameter', 'persistent_parameter'),
    )

    param_cols = {
            'Salinity': ['CTDSAL', 'SALNTY'],
            'Oxygen': ['CTDOXY', 'OXYGEN'],
            'Temperature': ['CTDTMP', 'REFTMP']
    }
    cols = ['BTLNBR', 'CTDPRS'] + param_cols[st.session_state.selected_parameter]

    # set up the table
    config = {
            'BTLNBR': st.column_config.NumberColumn(
                    'Bottle', width=None, step=1, help='...of beer on the wall.',
            ),
            'CTDPRS': st.column_config.NumberColumn(
                    'Pressure', width=None, step=0.01, help="Can you handle it?",
            ),
            'CTDTMP': st.column_config.NumberColumn(
                    'CTD Temp', width=None, help="It's gettin' hot in here...",
            ),
            'REFTMP': st.column_config.NumberColumn(
                    'Ref Temp', width=None, help='One cool cucumber.',
            ),
            'CTDOXY': st.column_config.NumberColumn(
                    'CTD Oxy', width=None, help="Breathe it in.",
            ),
            'OXYGEN': st.column_config.NumberColumn(
                    'Bottle Oxy', width=None, help='Comin\' in the air tonight.',
            ),
            'CTDSAL': st.column_config.NumberColumn(
                    'CTD Salinity', width=None, help="Flaky sea salt.",
            ),
            'SALNTY': st.column_config.NumberColumn(
                    'Bottle Salinity', width=None, help='Want some lime and tequila with that?',
            ),
            'Residual': st.column_config.NumberColumn(
                    width=None, step=0.0001, help='What difference does it make?',
            ),
    }

    # caluclate residuals for display
    residuals = pd.Series(
            data[param_cols[st.session_state.selected_parameter][1]] -
            data[param_cols[st.session_state.selected_parameter][0]],
            name='Residual'
    )
    table_data = pd.concat([data[cols], residuals], axis=1).round(4)

    # draw the table
    st.dataframe(table_data, column_config=config, hide_index=True)

with page_right:
    # set up plot selection and draw the selection bar
    plot_options = ['Param vs Pressure', 'Residual Plot']
    if 'persistent_plot' not in st.session_state:
        st.session_state.selected_plot = plot_options[0]
    else:
        st.session_state.selected_plot = st.session_state.persistent_plot
    st.segmented_control(
            "Plot View",
            options=plot_options,
            selection_mode="single",
            key='selected_plot',
            on_change=store_persistent,
            args=('selected_plot', 'persistent_plot'),
    )

    # set up the param vs pressure plot

    # Bottle data
    # transform data...
    pvp_point_data = (
            table_data[
                ['CTDPRS', 'BTLNBR'] +
                param_cols[st.session_state.selected_parameter]
            ].melt(['CTDPRS', 'BTLNBR'], var_name='Parameter', value_name='measured')
    )

    pvp_points = alt.Chart(pvp_point_data).mark_point(filled=True).encode(
            x=alt.X(
                    'measured:Q',
                    title=st.session_state.selected_parameter,
                    scale=alt.Scale(zero=False),
                    axis=alt.Axis(tickCount=9),
            ),
            y=alt.Y(
                    'CTDPRS:Q',
                    title='Pressure',
                    scale=alt.Scale(reverse=True),
                    axis=alt.Axis(tickCount=6),
            ),
            color=alt.Color('Parameter', scale=alt.Scale(range=['blue', 'red'])),
            shape=alt.Shape('Parameter', scale=alt.Scale(range=['circle', 'cross'])),
            size=alt.Size('Parameter', scale=alt.Scale(range=[100, 180])),
            opacity=alt.Opacity('Parameter', scale=alt.Scale(range=[0.8, 0.5])),
            tooltip=[
                    alt.Tooltip(field='BTLNBR', title='Bottle'),
                    alt.Tooltip(field='measured', title=st.session_state.selected_parameter),
                    alt.Tooltip(field='Parameter'),
            ]
    ).properties(
            width=500,
            height=700,
    ).interactive(
    )
    # continuous data...
    cast_data = pd.read_csv(
            Path(datadir, '%s.csv' % file_list[cast_list.index(st.session_state.selected_cast)]),
            header=0, skipfooter=1,
            skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],
            engine='python',
    )
    line_data = cast_data[
        ['CTDPRS', param_cols[st.session_state.selected_parameter][0]]
    ].melt('CTDPRS', var_name='Parameter', value_name='measured')

    pvp_lines = alt.Chart(line_data).mark_line(color='darkgrey', strokeWidth=1).encode(
            x=alt.X('measured:Q'),
            y=alt.Y('CTDPRS:Q'),
            order='CTDPRS',
    )

    param_vs_p = alt.layer(pvp_points, pvp_lines).configure_legend(
            title=None,
            strokeColor='gray',
            fillColor='#EEEEEE',
            padding=10,
            cornerRadius=4,
            orient='bottom-right'
    )

    # set up residual vs pressure plot
    resid_point_data = table_data[['CTDPRS', 'BTLNBR', 'Residual']]

    resid_vs_p = alt.Chart(resid_point_data).mark_point(size=50, filled=True).encode(
            x=alt.X(
                    'Residual:Q',
                    title='%s Residual' % st.session_state.selected_parameter,
                    scale=alt.Scale(
                            domain=st.session_state.range_limits[st.session_state.selected_parameter],
                            zero=False
                    ),
                    # axis=alt.Axis(tickCount=9),
            ),
            y=alt.Y(
                    'CTDPRS:Q',
                    title='Pressure',
                    scale=alt.Scale(reverse=True),
                    axis=alt.Axis(tickCount=6),
            ),
            color=alt.Color('Residual').scale(
                    domain=st.session_state.range_limits[st.session_state.selected_parameter],
                    range=['darkred', 'orangered', 'teal', 'orangered', 'darkred']
            ),
            tooltip=[
                    alt.Tooltip(field='BTLNBR', title='Bottle'),
                    alt.Tooltip(field='Residual'),
            ]
    ).properties(
            width=500,
            height=700,
    ).interactive()

    resid_limits = pd.DataFrame.from_dict(st.session_state.resid_limits[st.session_state.selected_parameter])

    residual_bounds_lines = pd.concat([
            resid_limits,
            pd.DataFrame({
                    'x': resid_limits['x'] * -1,
                    'y1': resid_limits['y1'],
                    'y2': resid_limits['y2']
            })
    ])

    resid_bounds = alt.Chart(residual_bounds_lines).mark_rule(
            color='darkgrey'
    ).encode(
            x='x',
            y='y1',
            y2='y2',
    )

    if st.session_state.selected_plot == 'Param vs Pressure':
        st.altair_chart(param_vs_p, use_container_width=False)
    elif st.session_state.selected_plot == 'Residual Plot':
        st.altair_chart(resid_vs_p + resid_bounds, use_container_width=False)
