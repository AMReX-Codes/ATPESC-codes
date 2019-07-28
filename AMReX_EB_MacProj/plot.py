#!/usr/bin/env python3
import yt
from yt.frontends.boxlib.data_structures import AMReXDataset
import argparse
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infiles', type=str, default="\Aplt[0-9]*\Z", help='Regex for input files to plot.')
args = parser.parse_args()

def doit(filename):
    ds = AMReXDataset(filename)
    s = yt.SlicePlot(ds, 'z', 'vfrac',origin="native")
    s.set_log("vfrac", 0)
    s.annotate_particles(1.0, p_size = 10.0)
    s.save("{}.png".format(filename))

if __name__ == "__main__":
    regexp = re.compile(args.infiles)
    # Files in the current directory
    cwd_files = os.listdir()
    # Filter files by regexp and skip if {file}.skip exists
    subject_files = []
    for f in cwd_files:
        if regexp.match(f):
            subject_files.append(f)

    for f in subject_files:
        doit(f)
