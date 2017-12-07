"""
retrieve all data related to a particular frame
"""

import sys, os, shutil
import pickle, glob
import argparse

import util

base_files = ["fargo", "Jup.cfg", "this.par", "this.sh", "planet0.dat", "bigplanet0.dat"]
gas_files = ["gasdens%d.dat", "gasvrad%d.dat", "gasvtheta%d.dat"]
dust_files = ["gasddens%d.dat", "gasdvrad%d.dat", "gasdvtheta%d.dat"]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot dust density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')

    # Directory Selection
    parser.add_argument('--dir', dest = "directory", default = None,
                         help = 'source directory -- required (default: None)')

    # Dust Files
    parser.add_argument('--no_dust', dest = "dust", action = 'store_false', default = True,
                         help = 'copy dust files (default: yes)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

frame_range = util.get_frame_range(args.frames)
directory = args.directory
dust = args.dust

if directory is None:
	print "Error: --dir directory is required!"
	exit()

###############################################################################

def copy_file(src, file, dst = "."):
	"""copy file from src to dst"""

	# Trim trailing slash from directories
	if src[-1] == "/":
		src = src[:-1]
	if dst[-1] == "/":
		dst = dst[:-1]

	shutil.copyfile("%s/%s" % (src, file), "%s/%s" % (dst, file))

def copy_files(frame):
	for file in base_files:
		copy_file(directory, file)

	for file in gas_files:
		copy_file(directory, file % frame)

	if dust:
		for file in dust_files:
			copy_file(directory, file % frame)

###############################################################################

### Copy Files ###

for frame in frame_range:
	copy_files(frame)
