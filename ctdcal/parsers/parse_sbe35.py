#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parses SBE35 reference temperature data from Seaterm capture files.

Currently requires test casts or any other extraneous bottle fires to be
manually trimmed out of input files.

Capture files must be named with cast id, or start-end cast id for a range
of casts, with the default .cap file extension.

Can alternatively be called from the command line to parse files individually.
usage: parse_sbe35.py [-h] -cf CAST_FILE -sc CAST_NAME [-o OUTDIR] infile
"""
import argparse
import logging
import os
import re

from datetime import datetime as dt
from datetime import timedelta
from pathlib import Path

from ctdcal.common import validate_dir

# Constants
PATTERN = (  # Pattern to match a valid line
    r'\s*([0-9]+)' +                           # line number
    r'\s+([0-9]+\s[A-Z][a-z]{2}\s[0-9]{4}' +   # date
    r'\s+[0-9]{2}:[0-9]{2}:[0-9]{2})' +        # time
    r'\sbn\s=\s+([0-9]+)' +                    # bottle number
    r'\sdiff\s=\s+([0-9]+)' +                  # diff (unused)
    r'\sval\s=\s+([0-9]{6}\.[0-9]{1})' +       # val (unused)
    r'\st90\s=\s*([+-]?[0-9]+\.[0-9]+)'        # temperature (t90)
)
NEWLINE = r'(?:\r\n|\n)?'
TIMEDELTA = 30                 # Max time between bottles of same cast in minutes
LINELEN = 79


class Cast:
    """
    Parse and export individual casts from a raw sensor dump.
    """
    def __init__(self, cast_name):
        self.cast_name = cast_name
        self.sample_buffer = []
        self.last_bottle = None
        self.last_timestamp = None

    def detect_new(self, next_timestamp):
        if self.last_timestamp is None:
            return False
        elif (next_timestamp - self.last_timestamp) > timedelta(minutes=TIMEDELTA):
            return True
        else:
            return next_timestamp < self.last_timestamp

    def parse_cast(self, lines, i):
        for line in lines[i:]:
            next_timestamp = timestamp_from_line(line)
            if self.detect_new(next_timestamp):
                return i
            self.sample_buffer.append(line)
            self.last_timestamp = next_timestamp
            i += 1

    def export_cast_file(self, outdir):
        """
        Write the sample buffer to a text file and name it according to the
        cast name.
        """
        fname = Path(outdir, '%s.csv' % self.cast_name)
        with open(fname, 'w') as outfile:
            for line in self.sample_buffer:
                outfile.write(','.join(line))
                outfile.write('\n')
        return True


def timestamp_from_line(line):
    # line = line.split()
    return dt.strptime(line[1], '%d %b %Y %H:%M:%S')


def parse_lines(lines, pattern):
    # Read input lines and parse out samples and cast data.
    parsed = []
    for line in lines:
        if len(line) != LINELEN:
            # skip lines with invalid length
            continue
        match = re.match(pattern, line)
        if match:
            parsed.append(match.groups())
        else:
            logging.warning("Line does not match the pattern.")
            continue
    return parsed


def parse_sbe35(casts, indir, outdir):
    """
    Searches indir for SBE35 data capture files, parses valid data lines by cast and
    exports as csv.
    
    Parameters
    ----------
    casts : list[str]
        sequential cast names
    indir : str or PathLike
        input directory
    outdir : str or PathLike
        output directory
    """
    # validate required directories
    indir = validate_dir(indir)
    outdir = validate_dir(outdir, create=True)

    # find raw input files
    raw_files = tuple(Path(indir).glob("*.cap"))
    infile_list = [f.stem for f in raw_files]
    if len(infile_list) < 1:
        raise FileNotFoundError('No SBE35 raw data files found in %s.' % indir)

    for i, infile in enumerate(raw_files):
        # parse filenames for starting and ending casts
        bounds = infile_list[i].split('-')
        starting_cast = casts.index(bounds[0])
        ending_cast = casts.index(bounds[-1])

        infile_casts = casts[starting_cast:ending_cast + 1]

        with open(infile, 'rb') as f:
            lines = [line.decode(errors='ignore') for line in f.readlines()]
        pattern = re.compile(PATTERN + NEWLINE, re.DOTALL)
        lines = parse_lines(lines, pattern)

        # Parse the casts and export as individual files...
        i = 0
        for cast_name in infile_casts:
            cast = Cast(cast_name)
            i = cast.parse_cast(lines, i)
            cast.export_cast_file(outdir)
            if i is None:
                break


def main():
    parser = argparse.ArgumentParser(description="""Splits the log file from
                                     an SBE35 into individual cast files.""")
    parser.add_argument('-cf', '--castfile', dest='cast_file', type=str,
                        help="""file containing the list of cast ids
                                (ssscc file)""", required=True)
    parser.add_argument('-sc', '--starting', dest='cast_name', type=str,
                        help='starting cast name', required=True)
    parser.add_argument('-o', '--outdir', dest='outdir', type=str,
                        default='split',
                        help='directory to save output files (default: ./split)')
    parser.add_argument('infile', type=str,
                        help='source file')

    args = parser.parse_args()

    infile = os.path.abspath(args.infile)
    cast_file = os.path.abspath(args.cast_file)
    outdir = os.path.abspath(args.outdir)
    cast_name = args.cast_name

    # parse cast file into a list of cast ids...
    with open(cast_file, 'r') as f:
        lines = f.readlines()
        casts = [c.strip() for c in lines]
        casts = casts[casts.index(cast_name):]

    # Test for output directory...
    validate_dir(outdir, create=True)

    # Read in input file. Sometimes serial line noise is recorded as
    # non-ascii characters in the file, so we'll open the file in binary
    # mode just in case...
    with open(infile, 'rb') as f:
        lines = [line.decode(errors='ignore') for line in f.readlines()]
    pattern = re.compile(PATTERN + NEWLINE, re.DOTALL)
    lines = parse_lines(lines, pattern)

    # Parse the casts and export as individual files...
    i = 0
    for cast_name in casts:
        cast = Cast(cast_name)
        i = cast.parse_cast(lines, i)
        cast.export_cast_file(outdir)
        if i is None:
            break


if __name__ == '__main__':
    main()
