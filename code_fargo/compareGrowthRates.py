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
fontsize = 17

my_dpi = 100
alpha = 0.5

max_frame = 0

# Set up figures
fig = plot.figure(figsize = (12, 6))
gs = gridspec.GridSpec(1, 2, width_ratios = [2, 1])
ax1 = plot.subplot(gs[0])
ax2 = plot.subplot(gs[1])

# Plot Curves
sq = 75
ax2.plot([-sq, sq], c = "k", linewidth = 2)

for taper in tapers:
    frame_range = pickle.load(open("excess_mass_frames_taper%d.p" % taper, "rb"))
    mass_over_time = pickle.load(open("excess_mass_values_taper%d.p" % taper, "rb"))

    # Log, Smooth, Take Derivative (to get growth rate)
    kernel_size = 1

    log_mass_over_time = np.log(mass_over_time) / np.log(10)
    smoothed_log_mass_over_time = smooth(log_mass_over_time, kernel_size)

    growth_rates = np.diff(smoothed_log_mass_over_time) / np.diff(frame_range)

    # Center growth rates on peak
    if taper == 10:
        max_index = 75
    else:
        max_index = np.argmax(growth_rates)
    #max_index = np.argmax(growth_rates)

    shifted_frames = np.array(frame_range) - frame_range[max_index]

    # Curves
    ax1.plot(frame_range, mass_over_time, linewidth = linewidth, label = r"$T_{growth}=$" + "%d" % taper)
    ax2.plot(shifted_frames[:-1], growth_rates, linewidth = linewidth)

    # Record Max Frame
    if frame_range[-1] > max_frame:
        max_frame = frame_range[-1]

# Annotate
title = r"$m_p = " + str(mass) + r" $ $M_J$, $\nu_{disk} = 10^{" + str(viscosity) + r"}$"
ax1.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
ax1.set_ylabel("Excess Mass", fontsize = fontsize)
ax1.set_title(title, fontsize = fontsize + 2)

ax2.set_xlabel(r"$t - t_{max-growth}$", fontsize = fontsize)
ax2.set_ylabel("Growth Rate", fontsize = fontsize)

ax1.legend(loc = "upper right", bbox_to_anchor = (2.1, 0.95))

# Axes
ax1.set_xlim(0.0, 1500)

ax2.set_xlim(-sq, sq)
ax2.set_xlim

# Save + Close
plot.savefig("growthRates_m%d_v%d.png" % (mass, abs(viscosity)))
plot.show()

plot.close()

