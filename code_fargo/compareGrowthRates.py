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
sq = 60 # Range on growth rate panel

for taper in tapers:
    frame_range = pickle.load(open("excess_mass_frames_taper%d.p" % taper, "rb"))
    mass_over_time = pickle.load(open("excess_mass_values_taper%d.p" % taper, "rb"))

    # Log, Smooth, Take Derivative (to get growth rate)
    kernel_size = 1

    log_mass_over_time = np.log(mass_over_time) / np.log(10)
    smoothed_log_mass_over_time = smooth(log_mass_over_time, kernel_size)

    growth_rates = np.diff(smoothed_log_mass_over_time) / np.diff(frame_range)

    # Center growth rates on peak
    #if taper == 10:
    #    max_index = np.searchsorted(frame_range, sq)
    #else:
    #    max_index = np.argmax(growth_rates)
    max_index = np.argmax(growth_rates)

    shifted_frames = np.array(frame_range) - frame_range[max_index]

    # Curves
    ax1.plot(frame_range, mass_over_time, linewidth = linewidth, label = r"$T_{growth}=$" + "%d" % taper)
    ax2.plot(shifted_frames[:-1], growth_rates, linewidth = linewidth)

    # Record Max Frame
    if frame_range[-1] > max_frame:
        max_frame = frame_range[-1]

# Reference Lines
ax1.plot([0, 10**(4)], [0.2, 0.2], c = "k", linewidth = 2)
ax2.plot([-sq, sq], [0, 0], c = "k", linewidth = 2)

# Annotate
title = r"$m_p = " + str(mass) + r" $ $M_J$, $\nu_{disk} = 10^{" + str(viscosity) + r"}$"
ax1.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
ax1.set_ylabel(r"$M_{excess}$", fontsize = fontsize)
ax1.set_title(title, y = 1.01, fontsize = fontsize + 2)

ax2.set_xlabel(r"$t - t_{max-growth}$", fontsize = fontsize)
ax2.set_ylabel(r"$dM_{excess}/dt$", fontsize = fontsize, labelpad = -8)
ax2.set_title("Growth Rates", fontsize = fontsize)

ax1.legend(loc = "lower right")

# Axes
ax1.set_xlim(0.0, 1500)
ax1.set_ylim(10**(-5), 4.0)
ax1.set_yscale("log")

ax2.set_xlim(-sq, sq)
ax2.set_ylim(-0.02, 0.10)

# For taper = 10 case only (if necessary)
#ax3 = ax2.twiny()
#ax3.set_xlim(0, 2 * sq)

# Save + Close
plot.savefig("growthRates_m%d_v%d.png" % (mass, abs(viscosity)), bbox_inches = "tight")
plot.savefig("growthRates_m%d_v%d.pdf" % (mass, abs(viscosity)), bbox_inches = "tight", format = "pdf")
plot.show()

plot.close()

