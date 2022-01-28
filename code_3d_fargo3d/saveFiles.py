"""
deletes specified files
"""

import sys, os, shutil
import glob, pickle
import argparse

#import numpy as np
#import util

###############################################################################

gas_files = ["gasdens%d.dat", "gasvx%d.dat", "gasvy%d.dat", "gasvz%d.dat", "gasenergy%d.dat"]
dust1_files = ["dust1dens%d.dat", "dust1vx%d.dat", "dust1vy%d.dat", "dust1vz%d.dat", "dust1energy%d.dat"]
dust2_files = ["dust2dens%d.dat", "dust2vx%d.dat", "dust2vy%d.dat", "dust2vz%d.dat", "dust2energy%d.dat"]
dust3_files = ["dust3dens%d.dat", "dust3vx%d.dat", "dust3vy%d.dat", "dust3vz%d.dat", "dust3energy%d.dat"]

gas_density_files = ["gasdens%d.dat"]
dust1_density_files = ["dust1dens%d.dat"]
dust2_density_files = ["dust2dens%d.dat"]
dust3_density_files = ["dust3dens%d.dat"]

gas_velocity_files = ["gasvx%d.dat", "gasvy%d.dat", "gasvz%d.dat"]
dust1_velocity_files = ["dust1vx%d.dat", "dust1vy%d.dat", "dust1vz%d.dat"]
dust2_velocity_files = ["dust2vx%d.dat", "dust2vy%d.dat", "dust2vz%d.dat"]
dust3_velocity_files = ["dust3vx%d.dat", "dust3vy%d.dat", "dust3vz%d.dat"]

density_files = ["gasdens%d.dat", "dust1dens%d.dat", "dust2dens%d.dat", "dust3dens%d.dat"]
velocity_files = ["gasvx%d.dat", "gasvy%d.dat", "gasvz%d.dat", "dust1vx%d.dat", "dust1vy%d.dat", "dust1vz%d.dat", "dust2vx%d.dat", "dust2vy%d.dat", "dust2vz%d.dat", "dust3vx%d.dat", "dust3vy%d.dat", "dust3vz%d.dat"]
energy_files = ["gasenergy%d.dat", "dust1energy%d.dat", "dust2energy%d.dat", "dust3energy%d.dat"]

all_files = []

def get_frame_range(frame_selection):
    """ return array of selected frames"""
    if len(frame_selection) == 1:
        frame_range = frame_selection
    elif len(frame_selection) == 2:
        start = frame_selection[0]; end = frame_selection[1]
        frame_range = range(start, end + 1)
    elif len(frame_selection) == 3:
        start = frame_selection[0]; end = frame_selection[1]; rate = frame_selection[2]
        frame_range = range(start, end + 1, rate)
    else:
        print "Error: Must supply 1, 2, or 3 frame arguments\nWith one argument, plots single frame\nWith two arguments, plots range(start, end + 1)\nWith three arguments, plots range(start, end + 1, rate)"
        exit()

    return frame_range

### Input Parameters ###

def new_argument_parser(description = "Plot convolved intensity maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')

    # File Selection
    parser.add_argument('--gas', dest = "gas", action = 'store_true', default = False,
                         help = 'delete gas-related (default: do not delete)')
    parser.add_argument('--dust1', dest = "dust1", action = 'store_true', default = False,
                         help = 'delete dust-related (default: do not delete)')
    parser.add_argument('--dust2', dest = "dust2", action = 'store_true', default = False,
                         help = 'delete dust-related (default: do not delete)')
    parser.add_argument('--dust3', dest = "dust3", action = 'store_true', default = False,
                         help = 'delete dust-related (default: do not delete)')

    parser.add_argument('--gasdens', dest = "gas_density", action = 'store_true', default = False,
                         help = 'delete gas density (default: do not delete)')
    parser.add_argument('--dust1dens', dest = "dust1_density", action = 'store_true', default = False,
                         help = 'delete dust density (default: do not delete)')
    parser.add_argument('--dust2dens', dest = "dust2_density", action = 'store_true', default = False,
                         help = 'delete dust density (default: do not delete)')
    parser.add_argument('--dust3dens', dest = "dust3_density", action = 'store_true', default = False,
                         help = 'delete dust density (default: do not delete)')

    parser.add_argument('--gasvxy', dest = "gas_velocity", action = 'store_true', default = False,
                         help = 'delete gas velocity (default: do not delete)')
    parser.add_argument('--dust1vxy', dest = "dust1_velocity", action = 'store_true', default = False,
                         help = 'delete dust velocity (default: do not delete)')
    parser.add_argument('--dust2vxy', dest = "dust2_velocity", action = 'store_true', default = False,
                         help = 'delete dust velocity (default: do not delete)')
    parser.add_argument('--dust3vxy', dest = "dust3_velocity", action = 'store_true', default = False,
                         help = 'delete dust velocity (default: do not delete)')

    parser.add_argument('--density', dest = "density", action = 'store_true', default = False,
                         help = 'delete density (default: do not delete)')
    parser.add_argument('--velocity', dest = "velocity", action = 'store_true', default = False,
                         help = 'delete velocity (default: do not delete)')
    parser.add_argument('--energy', dest = "energy", action = 'store_true', default = False,
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
save_range = get_frame_range(args.frames)
print "Save: (%d)" % len(save_range)
print save_range

all_files = range(args.frames[0], args.frames[1])
delete_range = [int(round(x, 0)) for x in all_files if x not in save_range]
print "Delete: (%d)" % len(delete_range)
print delete_range

###############################################################################

test_count = 0
delete_count = 0

save_count = 0

###############################################################################

def monitor_save(fns):
    """ count and list files to save """
    global save_count

    existing_fns = [fn for fn in fns if os.path.exists(fn)]
    print len(existing_fns) #, existing_fns

    save_count += len(existing_fns)
    print "Save Count: %d" % save_count

def trash(fns, delete):
    """ if delete, delete. else print. """
    global test_count
    global delete_count

    if delete:
        for fn in fns:
            if os.path.exists(fn):
                os.remove(fn)
                delete_count += 1
        print "Delete Count: %d" % delete_count
    else:
        existing_fns = [fn for fn in fns if os.path.exists(fn)]
        if len(existing_fns) > 1:
            print existing_fns[0], existing_fns[1], existing_fns[-1]
            print len(existing_fns) #, existing_fns
        elif len(existing_fns) == 0:
            print len(existing_fns), existing_fns
        
        test_count += len(existing_fns)
        print "Test Count: %d" % test_count

def gather_files_to_save(fn):
    """ gathers all such files in frame range """
    print fn
    fns = [fn % frame for frame in save_range]
    return fns

def gather_files_to_delete(fn):
    """ gathers all such files in frame range """
    print fn
    fns = [fn % frame for frame in delete_range]
    return fns

def delete_files():
    """ delete selected files """

    ###############################

    if args.gas:
        for fn in gas_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust1:
        for fn in dust1_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust2:
        for fn in dust1_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust3:
        for fn in dust1_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    ###############################

    if args.gas_density:
        for fn in gas_density_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust1_density:
        for fn in dust1_density_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust2_density:
        for fn in dust2_density_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust3_density:
        for fn in dust3_density_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    ###############################

    if args.gas_velocity:
        for fn in gas_velocity_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust1_velocity:
        for fn in dust1_velocity_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust2_velocity:
        for fn in dust2_velocity_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.dust3_velocity:
        for fn in dust3_velocity_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    ###############################

    if args.density:
        for fn in density_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.velocity:
        for fn in velocity_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    if args.energy:
        for fn in energy_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)

    ###############################

    if args.all:
        for fn in all_files:
            monitor_save(gather_files_to_save(fn))
            trash(gather_files_to_delete(fn), args.delete)


###############################################################################

#### Delete! ####

delete_files()


