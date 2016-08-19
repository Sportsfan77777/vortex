"""
compares excess masses with different taper times

Usage:
python compareExcessMassesOverTime.py
"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool
from multiprocessing import Array as mp_array

import math
import numpy as np
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot
from matplotlib import gridspec

from itertools import groupby

from pylab import rcParams # replace with rc ???
from pylab import fromfile

import util
from readTitle import readTitle

mass = 1
viscosity = -7
tapers = [10, 500, 1000, 2000]

# Helper Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

##### PLOTTING #####

# Plot Parameters
linewidth = 4
fontsize = 14

my_dpi = 100
alpha = 0.5

max_frame = 0

# Plot Curves
for taper in tapers:
    frame_range = pickle.load(open("excess_mass_frames_taper%d.p" % taper, "rb"))
    mass_over_time = pickle.load(open("excess_mass_values_taper%d.p" % taper, "rb"))

    # Log, Smooth, Take Derivative (to get growth rate)
    kernel_size = 1

    log_mass_over_time = np.log(mass_over_time) / np.log(10)
    smoothed_log_mass_over_time = smooth(log_mass_over_time, kernel_size)

    gradient_over_time = np.diff(smoothed_log_mass_over_time) / np.diff(frame_range)

    # Curve
    plot.plot(frame_range[:-1], gradient_over_time, linewidth = linewidth, label = r"$T_{growth}=$" + "%d" % taper)

    # Record Max Frame
    if frame_range[-1] > max_frame:
        max_frame = frame_range[-1]

# Annotate
title = r"$m_p = " + str(mass) + r" $ $M_J$, $\nu_{disk} = 10^{" + str(viscosity) + r"}$"
plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
plot.ylabel("Growth Rate", fontsize = fontsize)
plot.title(title, fontsize = fontsize + 2)

plot.legend(loc = "upper right")

# Axes
plot.xlim(0.0, 1500)
#plot.ylim(0, 2)

#plot.xscale("log")
#plot.yscale("log")

# Save + Close
plot.savefig("gradientExcessMassesOverTime_kernel%d_m%d_v%d.png" % (kernel_size, mass, abs(viscosity)))
plot.show()

plot.close()

