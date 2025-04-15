#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add and edit manual flags in a web browser.
"""
import pandas as pd
import streamlit as st

from ctdcal.fitting.fit_common import df_node_to_BottleFlags, save_node, bf_from_json
from ctdcal.tools.qc_demo.qc_demo_app import store_persistent


# begin drawing the page...
st.set_page_config(page_title='Edit Flags', page_icon=':ocean:', layout='wide')
st.subheader("Edit Manual Flags")


def update_rows():
    node = df_node_to_BottleFlags(
            st.session_state.edited_data.sort_values(
                    by=['cast_id', 'bottle_num'],
            ).astype(
                    {'cast_id': str, 'bottle_num': int, 'value': int, 'notes': str}
            )
    )
    save_node(st.session_state.flag_file, node, st.session_state.selected_flags, create_new=True)
    st.session_state.flag_data = bf_from_json(st.session_state.flag_file)


# load from manual flag file
if 'flag_data' not in st.session_state:
    # using session state prevents input loss from reloading the json before the save
    st.session_state.flag_data = bf_from_json(st.session_state.flag_file)
    st.session_state.edited_data = st.session_state.flag_data

# column configs for the flag node tables
config = {
        'cast_id': st.column_config.TextColumn(
                'Cast', width=None, help='A spell, maybe?', required=True
        ),
        'bottle_num': st.column_config.NumberColumn(
                'Bottle', width=None, help="Mine's full of happiness.", required=True
        ),
        'value': st.column_config.NumberColumn(
                'Flag', width=None, help="Capture it!", required=True
        ),
        'notes': st.column_config.TextColumn(
                'Comments', width='medium', help='Can you describe the ruckus?', required=True
        ),
}

parameter_list = [name for name in st.session_state.flag_data]
if 'persistent_flags' not in st.session_state:
    st.session_state.selected_flags = parameter_list[0]
else:
    st.session_state.selected_flags = st.session_state.persistent_flags

st.segmented_control(
        "Parameter",
        options=parameter_list,
        selection_mode="single",
        key='selected_flags',
        on_change=store_persistent,
        args=('selected_flags', 'persistent_flags'),
)

# draw a table for the selected flag node and make it editable...
st.markdown("**%s flags**" % st.session_state.selected_flags.capitalize())
flags = pd.DataFrame.from_dict(
        st.session_state.flag_data[st.session_state.selected_flags]).astype(
        {'cast_id': str, 'bottle_num': int, 'value': int, 'notes': str}
)

st.session_state.display_data['added_rows'] = [{"cast_id":"09909","bottle_num":3,"value":9,"notes":"Nine"}]
st.session_state.edited_data = st.data_editor(
        flags,
        key='display_data',
        num_rows="dynamic",
        column_order=['cast_id', 'bottle_num', 'value', 'notes'],
        column_config=config,
        width=550,
        hide_index=True,
)

# save the edits
update_rows()

# draw the form to add a new node (new parameter to enter flags for)...
col1, col2 = st.columns([1, 2])
# Use 1:2 (width ratio) columns so the form isn't as wide as the entire page
col1.markdown("**Add a new parameter:**")
# simple form with one text box and submit button
with col1.form('new_node', clear_on_submit=True, enter_to_submit=False):
    new_node = st.text_input('Parameter name', placeholder='None')
    submit = st.form_submit_button('Add')

# if user has filled out the box and clicks submit, add the new node and
# reload the page to see it
if submit and new_node != '':
    keys = ['cast_id', 'bottle_num', 'value', 'notes']
    st.session_state.flag_data.add_node(new_node, keys)
    st.session_state.flag_data.save(st.session_state.flag_file)
    st.rerun()

st.write(st.session_state.display_data)
