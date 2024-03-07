#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import re

from datetime import datetime as dt
from datetime import timedelta


class Cast:
    """Simple class to store and process cast details. Must be initialized
    with a list of cast ids (names).
    """
    def __init__(self, casts, outdir):
        self.casts = casts
        self.cast_id = casts[0]
        self.outdir = outdir
        self.sample_buffer = []

    def new_cast(self):
        """Export to a file and reset the sample buffer, then get the next
        cast id.
        """
        if len(self.sample_buffer) != 0:
            self.export_cast_file()
            self.sample_buffer = []
            self.cast_id = self.casts[self.casts.index(self.cast_id)+1]
        return True

    def export_cast_file(self):
        """Write the sample buffer to a text file and name it according to the
        cast id.
        """
        fname = '%s/%s.cap' % (self.outdir, self.cast_id)  # ...even though I really dislike using .cap
        with open(fname, 'w') as outfile:
            for line in self.sample_buffer:
                outfile.write('%s' % line)
        return True


def main():
    parser = argparse.ArgumentParser(description="""Splits the log file from
                                     an SBE35 into individual cast files.""")
    parser.add_argument('-cf', '--castfile', dest='cast_file', type=str,
                        help="""file containing the list of cast ids
                                (ssscc file)""", required=True)
    parser.add_argument('-sc', '--starting', dest='cast_id', type=str,
                        help='starting cast id', required=True)
    parser.add_argument('-o', '--outdir', dest='outdir', type=str,
                        default='split',
                        help='directory to save output files (default: ./split)')
    parser.add_argument('infile', type=str,
                        help='source file')

    args = parser.parse_args()

    infile = os.path.abspath(args.infile)
    cast_file = os.path.abspath(args.cast_file)
    outdir = os.path.abspath(args.outdir)
    cast_id = args.cast_id

    # parse cast file into a list of cast ids...
    with open(cast_file, 'r') as f:
        lines = f.readlines()
        casts = [c.strip() for c in lines]
        casts = casts[casts.index(cast_id):]

    # initialize a cast object with the starting cast id...
    cast = Cast(casts, outdir)

    # Test for output directory...
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Pattern to match a valid line...
    pattern = r'bn.+?diff.+?t90'

    # Read input file and parse out samples and cast data. Sometimes serial
    # line noise is recorded as non-ascii characters in the file, so we'll
    # open the file in binary mode just in case...
    with open(infile, 'rb') as f:
        last_timestamp = dt(9999, 1, 1)
        # Decode from binary, ignoring the invalid characters, and test for
        # lines that are not samples...
        for line in f:
            line = line.decode(errors='ignore')
            if not re.search(pattern, line):
                continue
            # Check if we have the start of a new cast...
            words = line.split()
            timestamp = dt.strptime(' '.join(words[1:5]), '%d %b %Y %H:%M:%S')
            if words[7] == '1' and timestamp - last_timestamp > timedelta(minutes=15):
                cast.new_cast()
            cast.sample_buffer.append(line)
            last_timestamp = timestamp

    # Save the last cast and exit...
    if len(cast.sample_buffer) != 0:
        cast.export_cast_file()


if __name__ == '__main__':
    main()
