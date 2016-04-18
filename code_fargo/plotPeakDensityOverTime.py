"""
plots "weighted average of derivative around vortex" over time

Usage:
python plotVortexTrigger.py
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np
from scipy.signal import gaussian

from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
from readTitle import readTitle

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

### Helper Methods
def find_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start:])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    #print "Max", peak_rad, peak_density
    return peak_rad, peak_density

### Data ###
rate = 2
times = range(50, num_frames, rate)

# Planet Location
orbit_data = np.loadtxt(orbit_fn)
sm_axes = orbit_data[:, 2] # Planet Semi-Major Axis

# Radial Density Profiles
peaks = []
sigmas = []
zone_sigmas = []
for frame in times:
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    # Find Peak
    peak_rad, peak_density = find_peak(averagedDensity)

    # Zone off 10 scale heights around the vortex
    #spread = 5.0 * scale_height # half-width
    #start_index = np.searchsorted(rad, peak_rad - spread)
    #end_index = np.searchsorted(rad, peak_rad + spread)
    #density_in_vortex_zone = density[start_index : end_index, :]

    # Find Peak
    peak_in_vortex_zone = np.max(density_in_vortex_zone)

    # Find Sigmas
    #sigma = np.std(density[peak_rad, :])
    #zone_sigma = np.std(density_in_vortex_zone)

    # Track
    peaks.append(peak_in_vortex_zone)
    #sigmas.append(sigma)
    #zone_sigmas.append(zone_sigma)


##### PLOTTING #####

# Plot Parameters
alpha = 0.5
fontsize = 14
linewidth = 4

def make_plot():
    # Data
    xs = np.array(times)

    # Curves
    plot.plot(xs, peaks, c = "blue", linewidth = linewidth, label = r"$\rho_{peak}$")
    #plot.plot(xs, sigmas, c = "red", linewidth = linewidth - 1, alpha = alpha, label = r"$\sigma$")
    #plot.plot(xs, zone_sigmas, c = "green", linewidth = linewidth - 1, alpha = alpha, label = r"$\sigma_{zone}$")

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Peak Density", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    # Limits
    plot.xlim(xs[0], xs[-1])

    # Save + Close
    plot.savefig("peakDensityOverTime.png")
    plot.show()

    plot.close()


make_plot()