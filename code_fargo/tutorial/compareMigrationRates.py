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
directories = ["earth1", "earth4", "earth16", "jupiter2", "jupiter1", "saturn1", "saturn-half"]

###############################################################################

## Set file names ##
orbit_fn = "orbit0.dat"

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size, mode = 'nearest') # smoothing filter

# Plot Parameters
kernel_size = 20

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

fontsize = 14
linewidth = 3

start_time = 10 # start time
end_time = 95 # end time

def add_track(directory, index):
    fn = "../%s/%s" % (directory, orbit_fn)
    data = np.loadtxt(fn)

    times = data[:, 0] / (2 * np.pi) # Convert to orbital times
    sm_axes = data[:, 2] # Planet Semi-Major Axis

    dt = times[1] - times[0] # Note: output is constant
    smoothed_sm_axes = smooth(sm_axes, kernel_size)
    migration_rates = -(np.diff(smoothed_sm_axes) / dt) / smoothed_sm_axes[:-1] # -(da/dt) / a

    start = np.searchsorted(times, start_time)
    end = np.searchsorted(times, end_time)

    xs = times[start : end]; ys = migration_rates[start : end]
    plot.plot(xs, ys, linewidth = linewidth, label = directory, c = colors[index])

def make_plot():
    # Curves
    for i, directory in enumerate(directories):
        add_track(directory, i)

    # Annotate
    plot.title("Migration Rates", fontsize = fontsize + 2)
    plot.xlabel(r"$t$", fontsize = fontsize)
    plot.ylabel(r"$-\frac{1}{a} \frac{da}{dt}$", fontsize = fontsize)

    plot.legend(loc = "upper right")

    # Limits
    plot.xlim(0, 1.5 * end_time)
    #plot.ylim(min_y, max_y)

    # Save and Close
    plot.savefig("migrationRateComparison.png", bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot()
