#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add and edit manual flags in a browser.
"""
import argparse
from pathlib import Path

import pandas as pd
import streamlit as st

from ctdcal.fitting.fit_common import df_node_to_BottleFlags, save_node, bf_from_json


def main():
    # parse command line args
    parser = argparse.ArgumentParser(
            description='Add and edit manual flags in a browser'
    )
    parser.add_argument('infile', type=str,
                        help='path to the JSON flags file')
    args = parser.parse_args()
    if not args.infile:
        parser.print_help()
        return

    # load from manual flag file
    flag_file = Path(args.infile)
    if 'flag_data' not in st.session_state:
        # using session state prevents input loss from reloading the json before the save
        st.session_state.flag_data = bf_from_json(flag_file)

    # begin drawing the page...
    st.subheader("Edit Manual Flags")

    # column configs for the flag node tables
    config = {
            'cast_id': st.column_config.TextColumn(
                    'Cast', width=None, help='What?', required=True
            ),
            'bottle_num': st.column_config.NumberColumn(
                    'Bottle', width=None, help="It's the bottle number!", required=True
            ),
            'value': st.column_config.NumberColumn(
                    'Flag', width=None, help="The flag goes here!", required=True
            ),
            'notes': st.column_config.TextColumn(
                    'Comments', width=None, help='Spam!', required=True
            ),
    }

    # draw tables for each flag node and make them editable...
    for node_name in st.session_state.flag_data:
        st.markdown("**%s flags**" % node_name.capitalize())
        flags = pd.DataFrame.from_dict(
                st.session_state.flag_data[node_name]).astype(
                {'cast_id': str, 'bottle_num': int, 'value': int, 'notes': str}
        )
        display_data = st.data_editor(
                flags,
                num_rows="dynamic",
                column_config=config,
                width=550,
                hide_index=True
        )
        # save updates
        save_changes(display_data, node_name, flag_file)

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
        st.session_state.flag_data.save(flag_file)
        st.rerun()


def save_changes(node, node_name, fname):
    node_bf = df_node_to_BottleFlags(
            node.astype(
                    {'cast_id': str, 'bottle_num': int, 'value': int, 'notes': str}
            )
    )
    save_node(fname, node_bf, node_name, create_new=True)


if __name__ == "__main__":
    main()
