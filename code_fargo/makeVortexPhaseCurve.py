"""
plots vortex phase curve
uses slope to calculate vortex rotational period

Usage:
python makeVortexPhaseCurve.py
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

# Data
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
kernel_size = int(int(fargo_par["Nsec"]) / 10.0)

if (len(sys.argv) == 3):
    start = int(sys.argv[1])
    end = int(sys.argv[2])
else:
    start = 150
    end = 250
times = range(start, end)
vortex_phases = []
for i in times:
    # Get radial vortensity profile

    density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    vrad = (fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta))

    vorticity = curl(vrad, vtheta, rad, theta)
    vortensity = vorticity / normalized_density[1:, 1:]
    averaged_w = np.average(vortensity, axis = 1)

    # Find vortex

    outer_disk_start = np.searchsorted(rad, 1.2) # look for min vortensity beyond r = 1.2
    vortex_rad_outer_index = np.argmin(averaged_w[outer_disk_start:])

    vortex_rad_index = vortex_rad_outer_index + outer_disk_start
    vortex_theta_index = np.argmin(smooth(vortensity[vortex_rad_index, :], kernel_size))

    vortex_theta = theta[vortex_theta_index]
    if len(vortex_phases) > 0:
        previous_theta = vortex_phases[-1]
        while (previous_theta < vortex_theta):
            vortex_theta -= 2 * np.pi # should be less than previous theta
    vortex_phases.append(vortex_theta)

# Convert to degrees
vortex_phases = (180.0 / np.pi) * (np.array(vortex_phases))


##### PLOTTING #####

linewidth = 4
fontsize = 14

def make_plot():
    # Curves
    plot.plot(times, vortex_phases, linewidth = linewidth)

    # Annotate
    plot.title("Vortex Location", fontsize = fontsize + 2)
    plot.xlabel("Angle", fontsize = fontsize)
    plot.ylabel("Planet Orbit", fontsize = fontsize)

    # Limits
    plot.xlim(times[0], times[-1])

    # Save and Close
    plot.savefig("vortexPhaseCurve.png", bbox_inches = 'tight')
    plot.show()

    plot.cla()

gradient = np.gradient(vortex_phases)
print "Mean Slope: %.2f, Median Slope: %.2f" % (np.mean(gradient), np.median(gradient))

make_plot()