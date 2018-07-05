"""
calculate shifts relative to gas in 3 mm grains
default option is to also shift 3 mm grains away from the minimum
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

### Sizes ###
sizes = np.array([1.0, 0.3, 0.1, 0.3])
size_labels = ["cm", "hcm", "mm", "hmm"]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot azimuthal density profiles."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Shift Parameters
    parser.add_argument('-s', dest = "num_scale_heights", type = float, default = 4.0,
                         help = 'number of scale heights (default: 4.0)')
    parser.add_argument('--min', dest = "min_shift", type = float, default = 0.0,
                         help = 'minimum shift to check (default: 0)')
    parser.add_argument('--max', dest = "max_shift", type = float, default = 359.0,
                         help = 'maximum shift to check (default: 359)')
    parser.add_argument('--num', dest = "num_shifts", type = int, default = 360,
                         help = 'number of shifts to check (default: 360)')

    parser.add_argument('--ref', dest = "reference", type = int, default = "hcm",
                         help = 'reference density (default: hcm)')

    return parser

###############################################################################

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Shift Parameters
num_scale_heights = args.num_scale_heights
min_shift = args.min_shift
max_shift = args.max_shift
num_shifts = args.num_shifts

reference = args.reference

### Get Fargo Parameters ###
fargo_par = util.get_pickled_parameters(directory = "../%s-size" % reference)

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
surface_density_zero = fargo_par["Sigma0"]
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]
taper = fargo_par["MassTaper"]

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Helper Functions ###

def save_shifts(size, size_label, reference = "hcm"):
	density = util.read_data(directory = "../%s-size" % size_label)


# Make each call

for size_i, size_label_i in zip(sizes, size_labels):
	save_shifts(size_i, size_label_i, reference = reference)
