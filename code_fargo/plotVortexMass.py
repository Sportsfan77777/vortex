"""
adds up total mass in vortex (defined by vorticity between 0.0 and 0.2 and radius between 1.2 and 2.5)
in reality, the vortex should be at ~1.7 to 2.0
the vorticity thresholds are not perfect (density thresholds based on 50% of peak might work better)

Usage:
python plotVortexMass.py
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

## Set file names ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    orbit_fn = "orbit1.dat"
else:
    # fargo
    orbit_fn = "orbit0.dat"

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

surface_density = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

mass_taper = float(fargo_par["MassTaper"])

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

# Vortex Mass
def vortex_mass(radius, theta, density, vortensity):
    """ total mass contained in vortex """
    # Differentials
    d_rad = np.diff(radius)
    d_theta = np.diff(theta)

    d_rad = np.append(d_rad, d_rad[-1])
    d_theta = np.append(d_theta, d_theta[-1])

    area = radius[:, None] * np.outer(d_rad, d_theta)

    ### Zero out density outside of vortex ###
    # Do not include r < 1.2 or r > 2.5 regardless
    inner_disk_rad = 1.2
    inner_disk_index = np.searchsorted(radius, inner_disk_rad)

    outer_disk_rad = 2.5
    outer_disk_index = np.searchsorted(radius, outer_disk_rad)
    
    zoom_density = density[inner_disk_index : outer_disk_index]

    # Do not include vorticity < 0 or vorticity > 0.2 (non-vortex mostly falls into the latter)
    zoom_vortensity = vortensity[inner_disk_index : outer_disk_index]
    zoom_density[(zoom_vortensity < 0.0) || (zoom_vortensity > 0.2)] = 0

    # Multiply each density grid cell by its area
    zoom_area = area[inner_disk_index : outer_disk_index]
    vortex_mass_grid = zoom_density * zoom_area

    # Find total mass in vortex

    total_mass = np.sum(vortex_mass_grid)
    return total_mass

### Data ###
rate = 10
times = range(0, num_frames, rate)

vortex_masses = []
for frame in times:
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    vrad = (fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % frame).reshape(num_rad, num_theta))

    vorticity = curl(vrad, vtheta, rad, theta)
    vortensity = vorticity / density[1:, 1:]

    vortex_mass_i = vortex_mass(rad, theta, density, vortensity)
    vortex_masses.append(vortex_mass_i)


##### PLOTTING #####

# Plot Parameters
fontsize = 14
linewidth = 4

def make_plot():
    # Data
    xs = np.array(times)
    ys = np.array(vortex_masses)

    # Curves
    plot.plot(xs, ys, linewidth = linewidth)

    # Annotate
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Vortex Mass", fontsize = fontsize)
    plot.title("Mass Tapering: %d" % mass_taper, fontsize = fontsize + 1)

    # Limits
    plot.xlim(xs[0], xs[-1])

    # Save + Close
    plot.savefig("vortexMass.png")
    plot.show()

    plot.close()


make_plot()

