#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compare RBR vs SBE35 vs CTD
"""
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

# change these directories to point to your local directories then
# run: streamlit run rbr_comparison.py
rbrdir = '/Users/als026/data/i09n_2025/cruise_data/cnv/rbr'
sbedir = '/Users/als026/data/i09n_2025/cruise_data/cnv/reft'
ctddir = '/Users/als026/data/i09n_2025/cruise_data/btl/ctd'
#

cast_list = sorted([cast.stem[:5] for cast in Path(rbrdir).glob('*.csv')])

def store_persistent(kw, kp):
    st.session_state[kp] = st.session_state[kw]


# begin drawing the page
st.set_page_config(page_title='RBR Comparison', page_icon=':ocean:', layout='wide')
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

    # load the cast data
    rbrfile = Path(rbrdir, '%s_reft.csv' % st.session_state.selected_cast)
    sbefile = Path(sbedir, '%s_reft.csv' % st.session_state.selected_cast)
    ctdfile = Path(ctddir, '%s_btl_mean.pkl' % st.session_state.selected_cast)
    data = pd.read_pickle(ctdfile)[['btl_fire_num', 'CTDPRS', 'CTDTMP1']]
    try:
        data = pd.merge(data, pd.read_csv(rbrfile)[['btl_fire_num', 'REFTMP']], how='left', on='btl_fire_num').rename(columns={'REFTMP': 'rbr_tmp'})
    except FileNotFoundError:
        data['rbr_tmp'] = pd.NA
    try:
        data = pd.merge(data, pd.read_csv(sbefile)[['btl_fire_num', 'REFTMP']], how='left', on='btl_fire_num').rename(columns={'REFTMP': 'sbe_tmp'})
    except FileNotFoundError:
        data['sbe_tmp'] = pd.NA
    data.rename(columns={'btl_fire_num': 'btl'}, inplace=True)

    # make residual cols
    data['ctd-rbr'] = data['CTDTMP1'] - data['rbr_tmp']
    data['ctd-sbe'] = data['CTDTMP1'] - data['sbe_tmp']
    data['sbe-rbr'] = data['sbe_tmp'] - data['rbr_tmp']

    # set up the table
    config = {
            'btl': st.column_config.NumberColumn(
                    'Bottle', width=None, step=1,
            ),
            'CTDPRS': st.column_config.NumberColumn(
                    'Pressure', width=None, step=0.01,
            ),
            'CTDTMP1': st.column_config.NumberColumn(
                    'CTD Temp', width=None, step=0.0001,
            ),
            'rbr_tmp': st.column_config.NumberColumn(
                    'RBR Temp', width=None, step=0.0001,
            ),
            'sbe_tmp': st.column_config.NumberColumn(
                    'SBE35 Temp', width=None, step=0.0001,
            ),
            'ctd-rbr': st.column_config.NumberColumn(
                    'CTD-RBR', width=None, step=0.0001,
            ),
            'ctd-sbe': st.column_config.NumberColumn(
                    'CTD-SBE35', width=None, step=0.0001,
            ),
            'sbe-rbr': st.column_config.NumberColumn(
                    'SBE35-RBR', width=None, step=0.0001,
            ),
    }

    st.dataframe(data, hide_index=True, column_config=config)

# Right Column
with page_right:
    # set up plot selection and draw the selection bar
    plot_options = ['All Temp vs Pressure', 'CTD-RBR', 'CTD-SBE35', 'SBE35-RBR']
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
            data[
                ['CTDPRS', 'btl', 'rbr_tmp', 'sbe_tmp']
            ].melt(['CTDPRS', 'btl'], var_name='Parameter', value_name='measured')
    )

    pvp_points = alt.Chart(pvp_point_data).mark_point(filled=True).encode(
            x=alt.X(
                    'measured:Q',
                    # title=st.session_state.selected_parameter,
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
                    alt.Tooltip(field='btl', title='Bottle'),
                    alt.Tooltip(field='measured'),
                    # alt.Tooltip(field='measured', title=st.session_state.selected_parameter),
                    alt.Tooltip(field='Parameter'),
            ]
    ).properties(
            width=500,
            height=700,
    ).interactive(
    )
    # continuous data...
    line_data = data[
        ['CTDPRS', 'CTDTMP1']
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

    def residual_plot():
        # set up residual plots
        resid_params = {'CTD-RBR': 'ctd-rbr', 'CTD-SBE35': 'ctd-sbe', 'SBE35-RBR': 'sbe-rbr'}
        resid_label = resid_params[st.session_state.persistent_plot]
        resid_point_data = data[['CTDPRS', 'btl', resid_label]]

        resid_vs_p = alt.Chart(resid_point_data).mark_point(size=50, filled=True).encode(
                x=alt.X(
                        resid_label,
                        title='%s' % st.session_state.persistent_plot,
                        scale=alt.Scale(
                                domain=[-0.05, 0.05],
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
                color=alt.Color(resid_label).scale(
                        domain=[-0.05, 0.05],
                        range=['darkred', 'orangered', 'teal', 'orangered', 'darkred']
                ),
                tooltip=[
                        alt.Tooltip(field='btl', title='Bottle'),
                        alt.Tooltip(field=resid_label),
                ]
        ).properties(
                width=500,
                height=700,
        ).interactive()

        residual_bounds_dict = {
                'x': [0.02, 0.01, 0.005, 0.002],
                'y1': [0, 500, 1000, 2000],
                'y2': [500, 1000, 2000, 6000],
        }
        resid_limits = pd.DataFrame.from_dict(residual_bounds_dict)

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

        # draw the plot
        st.altair_chart(resid_vs_p + resid_bounds, use_container_width=False)

    # display the selected chart
    if st.session_state.selected_plot == 'All Temp vs Pressure':
        st.altair_chart(param_vs_p, use_container_width=False)
    elif st.session_state.selected_plot in ['CTD-RBR', 'CTD-SBE35', 'SBE35-RBR']:
        residual_plot()