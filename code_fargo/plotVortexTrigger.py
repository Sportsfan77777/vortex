"""
plots (peak density) / (pseudo-half-width) / (radial distance to planet)**2 over time

Usage:
python plotPeakStrength.py
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

    peak_rad = rad[outer_disk_start + peak_rad_outer_index]
    peak_density = averagedDensity[peak_rad]

    print "Max", peak_rad, peak_density
    return peak_rad, peak_density

def find_min(averagedDensity, peak_rad):
    try:
        outer_disk_start = np.searchsorted(rad, 1.0) # look for max radial density beyond r = 1.1
        outer_disk_end = np.searchsorted(rad, peak_rad)
        min_rad_outer_index = np.argmin(averagedDensity[outer_disk_start : outer_disk_end])

        min_rad = rad[outer_disk_start + min_rad_outer_index]
        min_density = averagedDensity[min_rad]

        print "Min", min_rad, min_density
        return min_rad, min_density
    except:
        # No Gap Yet
        return peak_rad, 0

def find_slope(averagedDensity, start_rad, end_rad):
    start = np.searchsorted(rad, start_rad)
    end = np.searchsorted(rad, end_rad)

    derivative_around_vortex = np.diff(averagedDensity[start : end + 1]) / np.diff(rad[start : end + 1])
    slope_magnitudes = np.abs(derivative_around_vortex)

    # Weight based on proximity to peak
    try:
        num_points = len(slope_magnitudes)
        sigma = num_points / 3
        gaussian_weights = gaussian(num_points, sigma)
        print "Gaussian", sigma, gaussian_weights

        mean_slope = np.average(slope_magnitudes, weights = gaussian_weights)
        print "Slope", mean_slope
        return mean_slope
    except:
        # No Gap Yet
        return 0

### Data ###
rate = 25
times = range(50, num_frames, rate)

# Planet Location
orbit_data = np.loadtxt(orbit_fn)
sm_axes = orbit_data[:, 2] # Planet Semi-Major Axis

# Radial Density Profiles
strengths = []
for frame in times:
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    averagedDensity = np.average(density, axis = 1) / surface_density

    # Find Peak
    peak_rad, peak_density = find_peak(averagedDensity)

    # Find Min
    min_rad, min_density = find_min(averagedDensity, peak_rad)

    # Get Amplitude
    amplitude = peak_density - min_density

    # Calculated Weighted Average of Derivative Between Min_Rad and (Peak_Rad + (Peak_Rad - Min_Rad))
    start_rad = min_rad
    end_rad = 2.0 * peak_rad - min_rad # Equal Width on Both Sides of Peak

    weighted_mean_derivative = find_slope(averagedDensity, start_rad end_rad)

    # Combine into metric
    strength_metric = amplitude * weighted_mean_derivative
    strengths.append(strength_metric)

    print frame, strength_metric
    print ""


##### PLOTTING #####

# Plot Parameters
fontsize = 14
linewidth = 4

def make_plot():
    # Data
    xs = np.array(times)
    ys = np.array(strengths)

    # Curves
    plot.plot(xs, ys, linewidth = linewidth)

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Vortex Trigger", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    # Limits
    plot.xlim(xs[0], xs[-1])

    # Save + Close
    plot.savefig("vortexTrigger.png")
    plot.show()

    plot.close()


make_plot()