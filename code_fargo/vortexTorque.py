"""
calculates torque contribution from vortex

Usage:
python vortexTorque.py
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np


import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util

## Check frame ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    frame = 0
else:
    # fargo
    frame = 1

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

max_frame = util.find_max_frame()
num_frames = max_frame + 1

# Curl function
def curl(v_rad, v_theta, rad, theta):
    """ z-component of the curl (because this is a 2-D simulation)"""
    ### Start Differentials ###
    d_rad = np.diff(rad)
    d_theta = np.diff(theta)

    dv_rad = np.diff(v_rad, axis = 1)
    dv_theta = np.diff(rad[:, None] * v_theta, axis = 0)
    ### End Differentials ###

    # z-Determinant
    partial_one = dv_theta / d_rad[:, None]
    partial_two = dv_rad / d_theta

    z_curl = (partial_one[:, 1:] - partial_two[1:, :]) / rad[1:, None]

    # Shift out of rotating frame (http://arxiv.org/pdf/astro-ph/0605237v2.pdf)
    if frame == 1:
        z_curl += 2

    return z_curl

#### Set Up and Binary Output (.npy) and Text Output (.dat) ####
## Binary ##
npy_fn = "vortexTorque.npy"
npy_file = open(npy_fn, 'wb')

binary_array = np.zeros((8, num_frames)) - 1 # initialize to -1
np.save(npy_file, binary_array)
npy_file.close()

## Text ##
dat_fn = "vortexTorque.dat"
dat_file = open(dat_fn, 'w')

# (1) Frame, (2) Total, (3) Inner, (4) Outer, 
# (5) InnerPositive, (6) InnerNegative, (7) OuterPositive, (8) OuterNegative
column_widths = 14 * np.ones(2, dtype = int)
column_widths[0] = 7

a = "Frame".center(column_widths[0])
b = "Phasee".center(column_widths[0])
c = "Vortex Torque".center(column_widths[1])
first_line = "%s %s %s\n" % (a, b, c)
dat_file.write(first_line)
dat_file.close()

# Data


#### Add vortex phase into here (make it a separate function --- put it in util)

for i in range(num_frames)
	density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    vrad = (fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta))

    vorticity = curl(vrad, vtheta, rad, theta)
    vortensity = vorticity / normalized_density[1:, 1:]

    # the vortex is any cell where the vortensity is between the thresholds of 0.0 and 0.2
    vortex_indices = np.where(vortensity < 0.2 and vortensity > 0.0)

    torqueMap = util.torque(rad, theta, density)
    vortexTorque = np.sum(torqueMap[vortex_indices])