"""
plots torque over time

 *** uses a smoothing function ***
"""

import sys
import os
import subprocess
import pickle

import numpy as np
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc

from scipy import signal as sig
from scipy.ndimage import filters as ff

import util

## Set file names ##
orbit_fn = "orbit0.dat"

### Get FARGO Parameters ###
fargo_par = util.get_pickled_parameters() # Retrieve parameters from *.par file

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt(orbit_fn)
select = range(0, len(data[:, -1]), rate)
times = (data[:, 0])[select] / (2 * np.pi) # Convert to orbital times
sm_axes = (data[:, 2])[select] # Planet Semi-Major Axis

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size, mode = 'nearest') # smoothing filter

# Plot Parameters
alpha = 0.2 # for non-smoothed curves
fontsize = 14
linewidth = 3

def make_plot():
    # Data
    kernel = 50.0
    xs = times; ys = sm_axes; smoothed_ys = smooth(ys, kernel)
    plot.plot(xs, ys, c = "blue", alpha = alpha, linewidth = linewidth)
    plot.plot(xs, smoothed_ys, c = "blue", linewidth = linewidth)

    # Annotate
    plot.title("Migration Track", fontsize = fontsize + 2)
    plot.xlabel(r"$t$", fontsize = fontsize)
    plot.ylabel(r"Planet Semimajor Axis ($a$)", fontsize = fontsize)

    #plot.legend(loc = "upper right")

    # Limits
    plot.xlim(0, xs[-1])
    plot.ylim(min(ys) - 0.002, max(ys) + 0.015)

    # Save and Close
    plot.savefig("migrationTrack.png", bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot()