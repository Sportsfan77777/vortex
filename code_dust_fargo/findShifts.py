"""
calculate shifts relative to gas in 3 mm grains
default option is to also shift 3 mm grains away from the minimum
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
from multiprocessing import Array as mp_array
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

    parser.add_argument('--ref', dest = "reference_label", default = "hcm",
                         help = 'reference density (default: hcm)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

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

reference_label = args.reference_label

### Get Fargo Parameters ###
fargo_par = util.get_pickled_parameters(directory = "../%s-size" % reference_label)

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

def get_shift(args):
    # Args
    i, frame, size_label, reference_label = args

    # Read Data
    density = util.read_gas_data(frame, fargo_par, directory = "../%s-size" % size_label)
    reference_density = util.read_gas_data(frame, fargo_par, directory = "../%s-size" % reference_label)

    shift, theta_shift = az.find_shift(density, reference_density, fargo_par, center = False, num_scale_heights = num_scale_heights, min_shift = min_shift, max_shift = max_shift, num_shifts = num_shifts)

    shift_array[i] = shift
    theta_shift_array[i] = theta_shift


def save_shifts(size_label, reference_label = "hcm"):
    # Collect Data
    if len(frame_range) == 1:
        get_shift((0, frame, size_label, reference_label))
    else:
        if num_cores > 1:
            pool_args = [(i, frame, size_label, reference_label) for (i, frame) in enumerate(frame_range)]

            p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
            p.map(get_shift, pool_args)
            p.terminate()
        else:
            for i, frame in enumerate(frame_range):
                get_shift((i, frame, size_label, reference_label))

    # Save Data
    shift_array_np = np.array(shift_array)
    theta_shift_array_np = np.array(theta_shift_array)

    pickle.dump(shift_array_np, open("../%s-size/shift_lookup.p" % size_label, 'w'))
    pickle.dump(theta_shift_array_np, open("../%s-size/theta_lookup.p" % size_label, 'w'))


### Make each call ###

for size_i, size_label_i in zip(sizes, size_labels):
    # Storage Arrays
    shift_array = mp_array("d", len(frame_range))
    theta_shift_array = mp_array("d", len(frame_range))

    # Save Data
    save_shifts(size_label_i, reference_label = reference_label)
