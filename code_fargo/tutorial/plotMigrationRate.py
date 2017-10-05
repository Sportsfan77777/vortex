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

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt(orbit_fn)
select = range(0, len(data[:, -1]), rate)
times = (data[:, 0])[select] / (2 * np.pi) # Convert to orbital times
sm_axes = (data[:, 2])[select] # Planet Semi-Major Axis

kernel_size = 50.0
dt = times[1] - times[0] # Note: output is constant
smoothed_sm_axes = smooth(sm_axes, kernel_size)
migration_rates = -np.diff(smoothed_sm_axes) / dt

# Plot Parameters
alpha = 0.2 # for non-smoothed curves
fontsize = 16
linewidth = 3

def make_plot():
    # Data
    xs = sm_axes[:-1]; ys = migration_rates
    plot.plot(xs, ys, c = "blue", linewidth = linewidth)

    # Analytic
    ys_a = np.max(migration_rates) * np.power(xs, -0.25)
    plot.plot(xs, ys_a, alpha = alpha, linestyle = "--", linewidth = linewidth - 1)

    # Annotate
    plot.title("Migration Rate", fontsize = fontsize + 2)
    plot.xlabel(r"$t$", fontsize = fontsize)
    plot.ylabel(r"$-\frac{da}{dt}$", fontsize = fontsize)

    #plot.legend(loc = "upper right")

    # Limits
    plot.xlim(0, xs[-1])
    plot.ylim(min(ys), max(ys))

    # Save and Close
    plot.savefig("migrationRate.png", bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot()