"""
deletes specified files
"""

import sys, os, shutil
import glob, pickle
import argparse

import numpy as np
import util

###############################################################################

gas_files = ["gasdens%d.dat", "gasvrad%d.dat", "gasvtheta%d.dat"]
dust_files = ["gasddens%d.dat", "gasdvrad%d.dat", "gasdvtheta%d.dat", "gasDiag1%d.dat", "gasDiag2%d.dat", "gasDiag3%d.dat", "gasTStop%d.dat"]

density_files = ["gasdens%d.dat, gasddens%d.dat"]
velocity_files = ["gasvrad%d.dat", "gasvtheta%d.dat", "gasdvrad%d.dat", "gasdvtheta%d.dat"]
diag_files = ["gasDiag1%d.dat", "gasDiag2%d.dat", "gasDiag3%d.dat", "gasTStop%d.dat"]

all_files = ["gasdens%d.dat", "gasvrad%d.dat", "gasvtheta%d.dat", "gasddens%d.dat", "gasdvrad%d.dat", "gasdvtheta%d.dat", "gasDiag1%d.dat", "gasDiag2%d.dat", "gasDiag3%d.dat", "gasTStop%d.dat"]

### Input Parameters ###

def new_argument_parser(description = "Plot convolved intensity maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')

    # File Selection
    parser.add_argument('--gas', dest = "gas", action = 'store_true', default = False,
                         help = 'delete gas-related (default: do not delete)')
    parser.add_argument('--dust', dest = "dust", action = 'store_true', default = False,
                         help = 'delete dust-related (default: do not delete)')

    parser.add_argument('--density', dest = "density", action = 'store_true', default = False,
                         help = 'delete density (default: do not delete)')
    parser.add_argument('--velocity', dest = "velocity", action = 'store_true', default = False,
                         help = 'delete velocity (default: do not delete)')
    parser.add_argument('--diag', dest = "diag", action = 'store_true', default = False,
                         help = 'delete diag and t_stop (default: do not delete)')

    parser.add_argument('--all', dest = "all", action = 'store_true', default = False,
                         help = 'delete all files (default: do not delete)')

    # Test out
    parser.add_argument('-d', dest = "delete", action = 'store_true', default = False,
                         help = 'delete, instead of just test out (default: do not delete)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

# Frames
frame_range = util.get_frame_range(args.frames)

###############################################################################

def trash(fns, delete):
    """ if delete, delete. else print. """
    if delete:
        for fn in fns:
            if os.path.exists(fn):
                os.remove(fn)
    else:
        existing_fns = []
        for fn in fns:
            if os.path.exists(fn):
                existing_fns += [fn]
        print existing_fns
        print len(existing_fns)

def gather_files(fn):
    """ gathers all such files in frame range """
    fns = []
    for i, frame in enumerate(frame_range):
        fns += [fn % frame]
    return fns

def delete_files():
    """ delete selected files """
    if args.gas:
        for fn in gas_files:
            trash(gather_files(fn), args.delete)

    if args.dust:
        for fn in dust_files:
            trash(gather_files(fn), args.delete)

    ###############################

    if args.density:
        for fn in density_files:
            trash(gather_files(fn), args.delete)

    if args.velocity:
        for fn in velocity_files:
            trash(gather_files(fn), args.delete)

    if args.diag:
        for fn in diag_files:
            trash(gather_files(fn), args.delete)

    ###############################

    if args.all:
        for fn in all_files:
            trash(gather_files(fn), args.delete)


###############################################################################

#### Delete! ####

delete_files()


