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

import colormaps as cmaps
plot.register_cmap(name = 'viridis', cmap = cmaps.viridis)

mass = 1
viscosity = -7
tapers = [10, 500, 1000, 2000]

# Helper Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

##### PLOTTING #####
# Colors
viridis = matplotlib.cm.get_cmap("viridis")

colors = {}
colors[10] = "grey"
colors[250] = viridis(0) #"firebrick"
colors[500] = viridis(0.25) #"gold"
colors[1000] = viridis(0.5) #"forestgreen"
colors[2000] = viridis(0.75) #"cornflowerblue"
colors[4000] = viridis(0.99) #"darkorchid"

# Plot Parameters
linewidth = 4
fontsize = 17

my_dpi = 100
alpha = 0.5

max_frame = 0

# Set up figures
fig = plot.figure(figsize = (12, 4))
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

    log_mass_over_time = np.log(mass_over_time)
    smoothed_log_mass_over_time = smooth(log_mass_over_time, kernel_size)

    growth_rates = np.diff(smoothed_log_mass_over_time) / (2 * np.pi * np.diff(frame_range))

    # Center growth rates on peak
    #if taper == 10:
    #    max_index = np.searchsorted(frame_range, sq)
    #else:
    #    max_index = np.argmax(growth_rates)
    test_growth_rates = growth_rates[:-1]
    if mass == 1 and viscosity == -7 and taper == 1000:
        test_growth_rates[:40] = 0
    if mass == 5 and viscosity == -7 and taper == 4000:
        test_growth_rates[:50] = 0
    max_index = np.argmax(test_growth_rates)

    shifted_frames = np.array(frame_range) - frame_range[max_index]

    # Curves
    ax1.plot(frame_range, mass_over_time, c = colors[taper], linewidth = linewidth, label = r"$T_\mathrm{growth}=$" + "%d" % taper, zorder = 5)
    ax1.scatter([frame_range[max_index + 1]], [mass_over_time[max_index + 1]], c = "k", marker = "s", s = 100, zorder = 50)
    ax2.plot(shifted_frames[:-1], growth_rates, c = colors[taper], linewidth = linewidth)

    # Record Max Frame
    if frame_range[-1] > max_frame:
        max_frame = frame_range[-1]

# Reference Lines
ax1.plot([0, 10**(4)], [0.2, 0.2], c = "k", linewidth = 2)
ax2.plot([-sq, sq], [0, 0], c = "k", linewidth = 2)

# Annotate
title = r"$M_p = " + str(mass) + r" $ $M_J$, $\nu = 10^{" + str(viscosity) + r"}$"
ax1.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
ax1.set_ylabel(r"$M_\mathrm{excess}$ $/$ $\Sigma_0 r_p^2$", fontsize = fontsize)
ax1.set_title(title, y = 1.01, fontsize = fontsize + 2)

ax2.set_xlabel(r"$t - t_\mathrm{max-growth}$", fontsize = fontsize)
ax2.set_ylabel(r"$d\log M_\mathrm{excess}/dt$", fontsize = fontsize, labelpad = -6)
ax2.set_title("Growth Rates", fontsize = fontsize)

ax1.legend(loc = "lower right")

# Axes
ax1.set_xlim(0.0, 1500)
ax1.set_ylim(10**(-5), 4.0)
ax1.set_yscale("log")

ax2.set_xlim(-sq, sq)
ax2.set_ylim(-0.01, 0.06)

# For taper = 10 case only (if necessary)
#ax3 = ax2.twiny()
#ax3.set_xlim(0, 2 * sq)

# Save + Close
plot.savefig("growthRates_m%d_v%d.png" % (mass, abs(viscosity)), bbox_inches = "tight")
plot.savefig("growthRates_m%d_v%d.pdf" % (mass, abs(viscosity)), bbox_inches = "tight", format = "pdf")
plot.show()

plot.close()

