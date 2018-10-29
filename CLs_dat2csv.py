#!/usr/bin/env python
""" CLs_dat2csv.py --- Define a CSV format for our CLs measurements :).
    Requirements:
    - Pandas
"""
from __future__ import (print_function, division)
from pandas import (read_csv, concat)
from numpy import (loadtxt)
from glob import (glob)
from sys import (argv)

if __name__ == '__main__':
    if len(argv) != 3:
        print("Usage:", argv[0], "<path to *.dat files (ending in /)>",
                "<output filename (CSV)>")
        exit(2)

    # Load data
    fn = glob("%scl_mock*.dat" % (argv[1]))
    if len(fn) < 1:
        print("Something went wrong :/")
        print("Does the passed path end in '/'?")
        exit(1)
    ell = loadtxt("%sell.dat" % argv[1])
    cl = [read_csv(f, sep=' ', index_col=None, header=None) for f in fn]
    for cell in cl:
        cell.columns = ['zbin%d' % (iz + 1)
                for iz in range(cell.shape[1])]
    
    # construct structured DataFrame
    keys = ['mock%d' % (im) for im in range(len(cl))]
    df = concat(cl, keys=keys, axis=1)
    df.index = ell
    df.index.name = 'elleff'
    df.columns.names = ['data', 'zbin']
    df.swaplevel(i=0, j=1, axis=1).to_csv(argv[2])
