#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Entry point for "QC Demo App"
"""
from pathlib import Path

import streamlit as st
from munch import munchify

st.session_state.spam = 'Lovely spam!'
st.session_state.exchange_dir = Path('/Users/als026/data/i09n_2025/cruise_data/save')
st.session_state.flag_file = Path('/Users/als026/data/i09n_2025/cruise_data/flag/bottle_flags_manual.json')
st.session_state.range_limits = {
        'Salinity': [-0.03, 0.03],
        'Oxygen': [-10.0, 10.0],
        'Temperature': [-0.05, 0.05],
}
st.session_state.resid_limits = munchify({
        'Salinity': {
                'x': [0.02, 0.01, 0.005, 0.002],
                'y1': [0, 500, 1000, 2000],
                'y2': [500, 1000, 2000, 6000],
        },
        'Oxygen': {
                'x': [7.5, 5.0, 2.5, 1.5],
                'y1': [0, 500, 1000, 2000],
                'y2': [500, 1000, 2000, 6000],
        },
        'Temperature': {
                'x': [0.02, 0.01, 0.005, 0.002],
                'y1': [0, 500, 1000, 2000],
                'y2': [500, 1000, 2000, 6000],
        },
})


# functions
def store_persistent(kw, kp):
    st.session_state[kp] = st.session_state[kw]


pg = st.navigation(
        [st.Page('qc_plot.py', title='QC Plot'),
         st.Page('edit_flags.py', title='Edit Flags')]
)
pg.run()
