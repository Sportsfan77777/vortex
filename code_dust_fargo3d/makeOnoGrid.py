"""
plot azimuthally averaged density
then, makes movies

Usage:
python plotAveragedDensity.py frame_number <== plot one frame
python plotAveragedDensity.py -m <== make a movie instead of plots (plots already exist)
python plotAveragedDensity.py -1 <<<===== Plots a sample
python plotAveragedDensity.py <== plot all frames and make a movie

"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np

import matplotlib
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

#import util
#import azimuthal as az
from readTitle import readTitle

from advanced import Parameters
from reader import Fields

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

#### Subtracting the Initial Profile ####

# Initial Vortices

initial_widths_h4 = [0.034, 0.037, 0.037, 0.039]
initial_amplitudes_h4 = [0.182, 0.198, 0.217, 0.300]

initial_widths_h6 = [0.052, 0.058, 0.057, 0.063]
initial_amplitudes_h6 = [0.139, 0.167, 0.206, 0.201]

initial_widths_h8 = [0.081, 0.077, 0.086, 0.091]
initial_amplitudes_h8 = [0.148, 0.195, 0.209, 0.236]

# Later Generation Vortices

second_widths_h4 = [0.041, 0.041, 0.040, 0.043, 0.044, 0.046, 0.046] # (A = 0.5, t = 795), (A = 0.167, t = 900)
second_amplitudes_h4 = [0.820, 0.779, 0.675, 0.691, 0.698, 0.689, 0.686] # then (A = 0.05: t = 940, t = 1235, t = 1520, t = 1800, t = 2130)

second_widths_h6 = [0.061, 0.064, 0.058, 0.061] # first two are from A = 0.05 (t = 2450, t = 3550)
second_amplitudes_h6 = [0.716, 0.667, 0.641, 0.641] # the rest are from A = 0.02 (t = 2470, t = 3080)

# Vortices that don't reform

death_widths_h4 = [0.034, 0.043, 0.044, 0.051]
death_amplitudes_h4 = [0.315, 0.681, 0.620, 0.593]

death_widths_h6 = [0.061, 0.064, 0.077]
death_amplitudes_h6 = [0.687, 0.687, 0.518]

#### Subtracting \Sigma = 0.667 ####

# Initial Vortices

initial_widths_h4 = []
initial_amplitudes_h4 = []

initial_widths_h6 = []
initial_amplitudes_h6 = []

initial_widths_h8 = []
initial_amplitudes_h8 = []

# Later Generation Vortices

second_widths_h4 = [0.041, 0.041, 0.040, 0.043, 0.044, 0.046, 0.046] # (A = 0.5, t = 795), (A = 0.167, t = 900)
second_amplitudes_h4 = [0.820, 0.779, 0.675, 0.691, 0.698, 0.689, 0.686] # then (A = 0.05: t = 940, t = 1235, t = 1520, t = 1800, t = 2130)

second_widths_h6 = [0.061, 0.064, 0.058, 0.061] # first two are from A = 0.05 (t = 2450, t = 3550)
second_amplitudes_h6 = [0.716, 0.667, 0.641, 0.641] # the rest are from A = 0.02 (t = 2470, t = 3080)

# Vortices that don't reform

death_widths_h4 = [0.034, 0.043, 0.044, 0.051]
death_amplitudes_h4 = [0.315, 0.681, 0.620, 0.593]

death_widths_h6 = [0.061, 0.064, 0.077]
death_amplitudes_h6 = [0.687, 0.687, 0.518]

###############################################################################


##### PLOTTING #####

linewidth = 3
labelsize = 17
fontsize = 19
markersize = 9
dpi = 100

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

colors = ["darkviolet", "royalblue", "sandybrown"]

def make_plot(show = False):
    # Figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Reference Lines
    def marginally_stable(x, h = 0.06):

        def option1(x):
            alpha = 325.0 * (h / 0.1)**(-3) 
            beta = 3.00
            return alpha * (x)**(beta)

        def option2(x):
            alpha = 1060000.0 * (h / 0.1)**(-3) 
            beta = 5.72
            return alpha * (x)**(beta)

        return np.where(x < 0.05, option1(x), option2(x))

    x_ref = np.logspace(-2, -1, 201)
    y_h4 = marginally_stable(x_ref, h = 0.04)
    y_h6 = marginally_stable(x_ref, h = 0.06)
    y_h8 = marginally_stable(x_ref, h = 0.08)
    ref_h4, = plot.plot(x_ref, y_h4, linewidth = linewidth, c = colors[0], label = r"$h = 0.04$")
    ref_h6, = plot.plot(x_ref, y_h6, linewidth = linewidth, c = colors[1], label = r"$h = 0.06$")
    ref_h8, = plot.plot(x_ref, y_h8, linewidth = linewidth, c = colors[2], label = r"$h = 0.08$")

    # Actual Values (Initial, Secondary, Stable)
    i_h4 = plot.scatter(initial_widths_h4, initial_amplitudes_h4, s = 100, c = colors[0], zorder = 100)
    i_h6 = plot.scatter(initial_widths_h6, initial_amplitudes_h6, s = 100, c = colors[1], zorder = 100)
    i_h8 = plot.scatter(initial_widths_h8, initial_amplitudes_h8, s = 100, c = colors[2], zorder = 100)

    s_h4 = plot.scatter(second_widths_h4, second_amplitudes_h4, s = 150, c = colors[0], zorder = 100, marker = "*", alpha = 0.5)
    s_h6 = plot.scatter(second_widths_h6, second_amplitudes_h6, s = 150, c = colors[1], zorder = 100, marker = "*", alpha = 0.5)
    #s_h8, = plot.scatter(second_widths_h8, second_amplitudes_h8, s = 150, c = colors[2], zorder = 100, marker = "")

    d_h4 = plot.scatter(death_widths_h4, death_amplitudes_h4, s = 100, c = colors[0], zorder = 100, marker = "D", alpha = 0.5)
    d_h6 = plot.scatter(death_widths_h6, death_amplitudes_h6, s = 100, c = colors[1], zorder = 100, marker = "D", alpha = 0.5)
    #d_h8, = plot.scatter(death_widths_h8, death_amplitudes_h8, s = 100, c = colors[2], zorder = 100, marker = "")

    # Axes
    plot.xlim(0.028, 0.1)
    plot.ylim(0.03, 3.5)
    plot.xscale("log"); plot.yscale("log")

    plot.xticks([0.04, 0.06, 0.08, 0.1], ["0.04", "0.06", "0.08", "0.10"])
    plot.yticks([0.1, 0.3, 1], ["0.1", "0.3", "1.0"])

    # Annotate
    plot.xlabel(r"$\Delta w$ $/$ $r$", fontsize = fontsize)
    plot.ylabel(r"$A$ $/$ $\Sigma_0$", fontsize = fontsize)
    plot.title(r"Stability of gap-edge bumps", y = 1.02, fontsize = fontsize + 2)

    plot.text(0.0295, 2.35, "Unstable", fontsize = fontsize)
    plot.text(0.078, 0.038, "Stable", fontsize = fontsize)
    
    first_legend = plot.legend([d_h6, s_h6, i_h6], ["RWI-stable", "Re-triggered Vortex", "Initial Vortex"], loc = "lower left", fontsize = fontsize - 2, scatterpoints = 1)
    second_legend = plot.legend([ref_h4, ref_h6, ref_h8], [r"$h = 0.04$", r"$h = 0.06$", r"$h = 0.08$"], loc = "upper right", fontsize = fontsize - 2)
    ax = plot.gca().add_artist(first_legend)


    #if version is None:
    #    save_fn = "%s_onoGrid.png" % ()
    #else:
    #    save_fn = "v%04d_%s_onoGrid.png" % (version)
    plot.savefig("onoGrid.png", bbox_inches = 'tight', dpi = dpi)
    plot.savefig("onoGrid.pdf", bbox_inches = 'tight', dpi = dpi, format = "pdf")

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = True)