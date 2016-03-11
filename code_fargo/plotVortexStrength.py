"""
plots (peak density) / (pseudo-half-width) / (radial distance to planet)**2 over time

Usage:
python plotVortexStrength.py
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

max_frame = util.find_max_frame()
num_frames = max_frame + 1

### Data ###
rate = 10
times = range(100, num_frames, rate)

# Planet Location
orbit_data = np.loadtxt(orbit_fn)
sm_axes = orbit_data[:, 2] # Planet Semi-Major Axis

# Radial Density Profiles
strengths = []
for frame in times:
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    averagedDensity = np.average(density, axis = 1) / surface_density

    # Find Peak
    outer_disk_start = np.searchsorted(rad, 1.2) # look for min vortensity beyond r = 1.2
    vortex_rad_outer_index = np.argmax(averagedDensity[outer_disk_start:])
    vortex_rad = rad[outer_disk_start + vortex_rad_outer_index]

    peak_density = averagedDensity[vortex_rad]

    # Get Pseudo-Half-Width (rather than find the half-width, take the avg value of the two points dr = 0.15 away)
    dr = 0.15
    inner_half_width = np.searchsorted(rad, vortex_rad - dr)
    outer_half_width = np.searchsorted(rad, vortex_rad + dr)
    avg_half_peak = 0.5 * (averagedDensity[inner_half_width] + averagedDensity[outer_half_width])

    # Get Distance from Vortex to Planet
    distance = vortex_rad - sm_axes[frame]

    # Combine into metric: (peak density) / (pseudo-half-width) / (radial distance to planet)**2
    strength_metric = peak_density / avg_half_peak / distance**2
    strengths.append(strength_metric)


##### PLOTTING #####

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 4


def make_plot():
	# Data
	xs = np.array(times)
	ys = np.array(strengths)

	# Curves
	plot.plot(xs, ys)

	# Annotate
	plot.xlim("Number of Planet Orbits")
	plot.ylabel("Vortex Peak Strength")

	# Limits
	plot.xlim(xs[0], xs[-1])

	# Save + Close
	plot.savefig("vortexPeakStrength.png")
	plot.show()

	plot.close()


make_plot()