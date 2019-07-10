#!/usr/bin/env python3
import yt
from yt.frontends.boxlib.data_structures import AMReXDataset
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input file to plot.')
parser.add_argument('-f', '--field', type=str, default='xvel', help='Name of field to plot.')
args = parser.parse_args()

if __name__ == "__main__":
    ds = AMReXDataset(args.infile)
    s = yt.SlicePlot(ds, 'z', args.field)
    s.annotate_particles(1.0)
    s.save("{}.png".format(args.infile))
