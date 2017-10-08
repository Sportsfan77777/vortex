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

## Choose directories ##
directories = ["earth1", "earth2", "earth4", "jupiter1", "jupiter2", "saturn1"]

###############################################################################

## Set file names ##
orbit_fn = "orbit0.dat"

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size, mode = 'nearest') # smoothing filter

# Plot Parameters
kernel_size = 20

fontsize = 14
linewidth = 3

def add_track(directory):
    fn = "../%s/%s" % (directory, orbit_fn)
    data = np.loadtxt(fn)

    times = data[:, 0] / (2 * np.pi) # Convert to orbital times
    sm_axes = data[:, 2] # Planet Semi-Major Axis

    xs = times; ys = smooth(sm_axes, kernel_size)
    plot.plot(xs, ys, linewidth = linewidth, label = directory)

def make_plot():
    # Curves
    for directory in directories:
        add_track(directory)

    # Annotate
    plot.title("Migration Track", fontsize = fontsize + 2)
    plot.xlabel(r"$t$", fontsize = fontsize)
    plot.ylabel(r"Planet Semimajor Axis ($a$)", fontsize = fontsize)

    plot.legend(loc = "upper right")

    # Limits
    #plot.xlim(0, xs[-1])
    #plot.ylim(min_y, max_y)

    # Save and Close
    plot.savefig("migrationTrackComparison.png", bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot()
