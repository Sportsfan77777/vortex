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

# Initial Vortices

initial_widths_h4 = [0.034, 0.037, 0.037, 0.039]
initial_amplitudes_h4 = [0.182, 0.198, 0.217, 0.300]

initial_widths_h6 = [0.052, 0.058, 0.057, 0.063]
initial_amplitudes_h6 = [0.139, 0.167, 0.206, 0.201]

initial_widths_h8 = [0.081, 0.077, 0.086, 0.091]
initial_amplitudes_h8 = [0.148, 0.195, 0.209, 0.236]

# Later Generation Vortices

second_widths_h4 = []
second_amplitudes_h4 = []

second_widths_h6 = [0.061, 0.064, 0.058, 0.061] # first two are from A = 0.05, the rest are from A = 0.02
second_amplitudes_h6 = [0.716, 0.667, 0.641, 0.641]

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
    plot.plot(x_ref, y_h4, linewidth = linewidth, c = colors[0])
    plot.plot(x_ref, y_h6, linewidth = linewidth, c = colors[1])
    plot.plot(x_ref, y_h8, linewidth = linewidth, c = colors[2])

    # Actual Values
    plot.scatter(initial_widths_h4, initial_amplitudes_h4, s = 100, c = colors[0], zorder = 100)
    plot.scatter(initial_widths_h6, initial_amplitudes_h6, s = 100, c = colors[1], zorder = 100)
    plot.scatter(initial_widths_h8, initial_amplitudes_h8, s = 100, c = colors[2], zorder = 100)

    #plot.scatter(second_widths_h4, second_amplitudes_h4, s = 150, c = colors[0], zorder = 100, marker = "*", alpha = 0.5)
    plot.scatter(second_widths_h6, second_amplitudes_h6, s = 150, c = colors[1], zorder = 100, marker = "*", alpha = 0.5)
    #plot.scatter(second_widths_h8, second_amplitudes_h8, s = 150, c = colors[2], zorder = 100, marker = "")

    plot.scatter(death_widths_h4, death_amplitudes_h4, s = 100, c = colors[0], zorder = 100, marker = "D", alpha = 0.5)
    plot.scatter(death_widths_h6, death_amplitudes_h6, s = 100, c = colors[1], zorder = 100, marker = "D", alpha = 0.5)
    #plot.scatter(death_widths_h8, death_amplitudes_h8, s = 100, c = colors[2], zorder = 100, marker = "")


    # Axes
    plot.xlim(0.028, 0.1)
    plot.ylim(0.008, 5)
    plot.xscale("log"); plot.yscale("log")

    # Annotate
    plot.xlabel(r"$\Delta w$ $/$ $r$", fontsize = fontsize)
    plot.ylabel(r"$A$ $/$ $\Sigma_0$", fontsize = fontsize)
    plot.title("RWI stability criteria", fontsize = fontsize + 2)

    #if version is None:
    #    save_fn = "%s_onoGrid.png" % ()
    #else:
    #    save_fn = "v%04d_%s_onoGrid.png" % (version)
    save_fn = "onoGrid.png" % ()
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = True)