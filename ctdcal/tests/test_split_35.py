#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:package: ctdcal.tests.test_split35
:file: ctdcal/tests/split_35.py
:author: Allen Smith
:brief: module unit tests
"""
from datetime import datetime as dt
from ctdcal.tools.split_35 import Cast, timestamp_from_line, parse_lines, PATTERN

CASTS = ['one', 'two', 'three']

test_data = ["0 2 Mar 2000  23:36:00 bn = 35",
             "1 2 Mar 2000  23:38:00 bn = 36",
             "2 2 Mar 2000  23:40:00 bn =  2",
             "3 3 Mar 2000  02:39:00 bn =  4"]
# test_data = [line.split() for line in raw]


def test_new_cast():
    # Tests four cases: first line of cast, sequential number same cast, non-sequential
    # number same cast, sequential number different cast.
    cast = Cast(CASTS[0])
    assert cast.detect_new(timestamp_from_line(test_data[0])) is False
    cast.last_timestamp = timestamp_from_line(test_data[0])
    assert cast.detect_new(timestamp_from_line(test_data[1])) is False
    assert cast.detect_new(timestamp_from_line(test_data[2])) is False
    assert cast.detect_new(timestamp_from_line(test_data[3])) is True


def test_parse_one_cast():
    cast = Cast(CASTS[0])
    cast.parse_cast(test_data, 0)
    assert len(cast.sample_buffer) == 3


def test_parse_all_cast():
    expected = [3, 1]
    buffer_lens = []
    i = 0
    for cast_name in CASTS:
        cast = Cast(cast_name)
        i = cast.parse_cast(test_data, i)
        buffer_lens.append(len(cast.sample_buffer))
        if i is None:
            break
    assert len(buffer_lens) == len(expected)
    assert buffer_lens[0] == expected[0]
    assert buffer_lens[1] == expected[1]


raw_input = """S>ds
SBE35 V 2.0a    SERIAL NO. 0001   01 Mar 2000  12:00:00
number of measurement cycles to average = 1
number of data points stored in memory = 3
bottle confirm interface = SBE 911plus
S>
S>
S>dd
  1 01 Mar 2024  05:00:00 bn =  1 diff =    11 val = 999.9 t90 =  9.999
  2 01 Mar 2024  05:01:00 bn =  2 diff =    11 val = 999.9 t90 =  9.999
  3 01 Mar 2024  05:02:00 bn =  3 diff =    11 val = 999.9 t90 =  9.999
S>
S>
"""


def test_parse_lines():
    lines = parse_lines(raw_input.splitlines(), PATTERN)
    assert len(lines) == 3
    assert timestamp_from_line(lines[0]) == dt.strptime('01 Mar 2024  05:00:00', '%d %b %Y %H:%M:%S')
