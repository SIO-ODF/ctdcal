#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:package: ctdcal.fitting.common
:file: ctdcal/fitting/common.py
:author: Allen Smith
:brief: Common code for use across ctdcal fitting modules
"""
import json
from pathlib import Path

from munch import Munch


class BottleFlags(Munch):
    """
    A dictionary class with the attribute-style access of Munch, plus methods
    for adding nodes and flag data, and loading or saving to/from a JSON file.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_node(self, label, keys):
        """
        Add a new empty node.

        Parameters
        ----------
        label (str) - name of node
        keys (list) - names of node keys
        """
        key_dict = {k: [] for k in keys}
        node = ({label: BottleFlags(key_dict)})
        self.update(node)

    def update_node(self, **kwargs):
        for k, v in kwargs.items():
            self[k].append(v)

    def save(self, fname):
        with open(fname, 'w') as f:
            f.write(self.toJSON())


# Function definitions
# --------------------

# BottleFlag wrangling
def df_node_to_BottleFlag(df, label):
    node_dict = df.to_dict()
    for k, v in node_dict.items():
        node_dict[k] = [vv for kk, vv in v.items()]
    return BottleFlags({label: BottleFlags(node_dict)})