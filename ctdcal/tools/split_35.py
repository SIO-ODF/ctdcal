#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import re

from datetime import datetime as dt
from datetime import timedelta

# Constants
PATTERN = r'bn.+?diff.+?t90'   # Pattern to match a valid line
TIMEDELTA = 15                 # Max time between bottles of same cast in minutes

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
        fname = '%s/%s.cap' % (outdir, self.cast_name)  # ...even though I really dislike using .cap
        with open(fname, 'w') as outfile:
            for line in self.sample_buffer:
                outfile.write('%s' % line)
        return True


def timestamp_from_line(line):
    line = line.split()
    return dt.strptime(' '.join(line[1:5]), '%d %b %Y %H:%M:%S')


def parse_lines(lines, pattern):
    # Read input lines and parse out samples and cast data.
    parsed = []
    for line in lines:
        if re.search(pattern, line):
            parsed.append(line)
    return parsed


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
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Read in input file. Sometimes serial line noise is recorded as
    # non-ascii characters in the file, so we'll open the file in binary
    # mode just in case...
    with open(infile, 'rb') as f:
        lines = [line.decode(errors='ignore') for line in f.readlines()]
    lines = parse_lines(lines, PATTERN)

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
