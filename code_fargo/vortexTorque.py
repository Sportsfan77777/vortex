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

from scipy.ndimage import filters as ff

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
num_columns = 3

## Binary ##
npy_fn = "vortexTorque.npy"
npy_file = open(npy_fn, 'wb')

binary_array = np.zeros((num_columns, num_frames)) - 1 # initialize to -1
np.save(npy_file, binary_array)
npy_file.close()

## Text ##
dat_fn = "vortexTorque.dat"
dat_file = open(dat_fn, 'w')

# (1) Frame, (2) Total, (3) Inner, (4) Outer, 
# (5) InnerPositive, (6) InnerNegative, (7) OuterPositive, (8) OuterNegative
column_widths = 14 * np.ones(num_columns, dtype = int)
column_widths[0] = 7
column_widths[1] = 10

a = "Frame".center(column_widths[0])
b = "Phase".center(column_widths[1])
c = "Vortex Torque".center(column_widths[2])
first_line = "%s %s %s\n" % (a, b, c)
dat_file.write(first_line)
dat_file.close()

# Data
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
kernel_size = int(int(fargo_par["Nsec"]) / 10.0)

for frame in range(num_frames):
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    vrad = (fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % frame).reshape(num_rad, num_theta))

    vorticity = curl(vrad, vtheta, rad, theta)
    vortensity = vorticity / normalized_density[1:, 1:]
    averaged_w = np.average(vortensity, axis = 1)

    ### Vortex Torque ###
    # the vortex is any cell where the vortensity is between the thresholds of 0.0 and 0.2
    #vortex_indices = np.where(vortensity < 0.2 and vortensity > 0.0)

    torqueMap = util.torque(rad, theta, density)
    vortexTorque = np.sum(torqueMap[(vortensity < 0.2) & (vortensity > 0.0)])

    ### Get Phase Too ###
    outer_disk_start = np.searchsorted(rad, 1.2) # look for min vortensity beyond r = 1.2
    vortex_rad_outer_index = np.argmin(averaged_w[outer_disk_start:])

    vortex_rad_index = vortex_rad_outer_index + outer_disk_start

    # Check 10 + 1 + 10 profiles
    dr = 10 # 10 indices
    vortex_width = range(vortex_rad_index - dr, vortex_rad_index + dr + 1)
    vortex_thetas = []
    for rad_index in vortex_width:
        vortex_theta_index = np.argmin(smooth(vortensity[rad_index, :], kernel_size))
        vortex_theta = theta[vortex_theta_index]
        vortex_thetas.append(vortex_theta)

    final_vortex_theta = np.median(vortex_thetas)

    # Format into strings
    scaling = 10**6 # multiply by one million to make things readable
    a = ("%d" % frame).center(column_widths[0])
    b = ("%.2f" % (final_vortex_theta)).center(column_widths[1])
    c = ("%.8f" % (vortexTorque * scaling)).center(column_widths[2])

    line = "%s %s %s\n" % (a, b, c)

    # Write to File
    dat_file = open(dat_fn, 'a')
    dat_file.write(line)
    dat_file.close()

    # Fill in entries
    binary_array[0, frame] = frame
    binary_array[1, frame] = final_vortex_theta
    binary_array[2, frame] = vortexTorque

    # Write to Binary
    npy_file = open(npy_fn, 'wb')
    np.save(npy_file, binary_array)
    npy_file.close()