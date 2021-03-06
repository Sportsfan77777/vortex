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

scale_height_version = True

###############################################################################

#### Subtracting the Initial Profile ####

# Initial Vortices

initial_accretion_numbers_h4 = np.array([1, 2, 3, 4])
initial_critical_maxima_h4 = np.array([1.120, 1.127, 1.130, 1.117])
initial_density_maxima_h4 = np.array([1.160, 1.160, 1.156, 1.146])

initial_accretion_numbers_h6 = np.array([1, 2, 3, 4])
initial_critical_maxima_h6 = np.array([1.168, 1.182, 1.200, 1.217])
initial_density_maxima_h6 = np.array([1.220, 1.228, 1.242, 1.257])

initial_accretion_numbers_h8 = np.array([1, 2, 3, 4])
initial_critical_maxima_h8 = np.array([1.261, 1.287, 1.316, 1.343])
initial_density_maxima_h8 = np.array([1.318, 1.348, 1.372, 1.389])

# Later Generation Vortices

second_accretion_numbers_h4 = np.array([2, 3, 3, 4, 4, 4])
second_critical_maxima_h4 = np.array([1.243, 1.192, 1.204, 1.172, 1.164, 1.174, 1.178])
second_density_maxima_h4 = np.array([1.343, 1.300, 1.310, 1.250, 1.257, 1.264, 1.275])

second_accretion_numbers_h6 = np.array([3, 3, 4, 4, 4])
second_critical_maxima_h6 = np.array([1.228, 1.228, 1.297, 1.250, 1.248])
second_density_maxima_h6 = np.array([1.396, 1.414, 1.354, 1.364, 1.375])

second_accretion_numbers_h8 = np.array([])
second_critical_maxima_h8 = np.array([])
second_density_maxima_h8 = np.array([])

# Vortices that reform interior

interior_accretion_numbers_h4 = np.array([4])
interior_critical_maxima_h4 = np.array([1.193])
interior_density_maxima_h4 = np.array([1.304])

interior_accretion_numbers_h6 = np.array([3, 3, 4])
interior_critical_maxima_h6 = np.array([1.224, 1.235, 1.304])
interior_density_maxima_h6 = np.array([1.432, 1.457, 1.382])

interior_accretion_numbers_h8 = np.array([])
interior_critical_maxima_h8 = np.array([])
interior_density_maxima_h8 = np.array([])

# Vortices that reform too interior

too_interior_accretion_numbers_h4 = np.array([1, 2, 3, 4])
too_interior_critical_maxima_h4 = np.array([1.202, 1.200, 1.189, 1.164])
too_interior_density_maxima_h4 = np.array([1.382, 1.357, 1.325, 1.304])

too_interior_accretion_numbers_h6 = np.array([1, 2])
too_interior_critical_maxima_h6 = np.array([1.264, 1.256])
too_interior_density_maxima_h6 = np.array([1.540, 1.465])

too_interior_accretion_numbers_h8 = np.array([])
too_interior_critical_maxima_h8 = np.array([])
too_interior_density_maxima_h8 = np.array([])

#### Subtracting \Sigma = 0.667 ####

# Initial Vortices

#initial_widths_h4 = []
#initial_amplitudes_h4 = []

#initial_widths_h6 = [0.061, 0.065, 0.066, 0.069]
#initial_amplitudes_h6 = [0.268, 0.313, 0.352, 0.348]

#initial_widths_h8 = []
#initial_amplitudes_h8 = []

# Later Generation Vortices

#econd_widths_h4 = [0.041, 0.041, 0.040, 0.043, 0.044, 0.046, 0.046] # (A = 0.5, t = 795), (A = 0.167, t = 900)
#second_amplitudes_h4 = [0.820, 0.779, 0.675, 0.691, 0.698, 0.689, 0.686] # then (A = 0.05: t = 940, t = 1235, t = 1520, t = 1800, t = 2130)

#second_widths_h6 = [0.061, 0.064, 0.058, 0.061] # first two are from A = 0.05 (t = 2450, t = 3550)
#second_amplitudes_h6 = [0.716, 0.667, 0.641, 0.641] # the rest are from A = 0.02 (t = 2470, t = 3080)

# Vortices that don't reform

#death_widths_h4 = [0.034, 0.043, 0.044, 0.051]
#death_amplitudes_h4 = [0.315, 0.681, 0.620, 0.593]

#death_widths_h6 = [0.061, 0.064, 0.077]
#death_amplitudes_h6 = [0.687, 0.687, 0.518]

if scale_height_version:
    initial_critical_maxima_h4 = (initial_critical_maxima_h4 - 1.0) / 0.04
    initial_critical_maxima_h6 = (initial_critical_maxima_h6 - 1.0) / 0.06
    initial_critical_maxima_h8 = (initial_critical_maxima_h8 - 1.0) / 0.08

    initial_density_maxima_h4 = (initial_density_maxima_h4 - 1.0) / 0.04
    initial_density_maxima_h6 = (initial_density_maxima_h6 - 1.0) / 0.06
    initial_density_maxima_h8 = (initial_density_maxima_h8 - 1.0) / 0.08

    second_critical_maxima_h4 = (second_critical_maxima_h4 - 1.0) / 0.04
    second_critical_maxima_h6 = (second_critical_maxima_h6 - 1.0) / 0.06
    second_critical_maxima_h8 = (second_critical_maxima_h8 - 1.0) / 0.08

    second_density_maxima_h4 = (second_density_maxima_h4 - 1.0) / 0.04
    second_density_maxima_h6 = (second_density_maxima_h6 - 1.0) / 0.06
    second_density_maxima_h8 = (second_density_maxima_h8 - 1.0) / 0.08

    interior_critical_maxima_h4 = (interior_critical_maxima_h4 - 1.0) / 0.04
    interior_critical_maxima_h6 = (interior_critical_maxima_h6 - 1.0) / 0.06
    interior_critical_maxima_h8 = (interior_critical_maxima_h8 - 1.0) / 0.08

    interior_density_maxima_h4 = (interior_density_maxima_h4 - 1.0) / 0.04
    interior_density_maxima_h6 = (interior_density_maxima_h6 - 1.0) / 0.06
    interior_density_maxima_h8 = (interior_density_maxima_h8 - 1.0) / 0.08

    too_interior_critical_maxima_h4 = (too_interior_critical_maxima_h4 - 1.0) / 0.04
    too_interior_critical_maxima_h6 = (too_interior_critical_maxima_h6 - 1.0) / 0.06
    too_interior_critical_maxima_h8 = (too_interior_critical_maxima_h8 - 1.0) / 0.08

    too_interior_density_maxima_h4 = (too_interior_density_maxima_h4 - 1.0) / 0.04
    too_interior_density_maxima_h6 = (too_interior_density_maxima_h6 - 1.0) / 0.06
    too_interior_density_maxima_h8 = (too_interior_density_maxima_h8 - 1.0) / 0.08

###############################################################################


##### PLOTTING #####

linewidth = 3
labelsize = 17
fontsize = 22
markersize = 9
alpha = 0.7
dpi = 100
size_offsets = np.array([75, 50, 25, 0])

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

    x_ref = np.linspace(0, 10, 201)
    y_ref = x_ref
    #ref, = plot.plot(x_ref, y_ref, linewidth = linewidth, c = 'k')

    if scale_height_version:
        # Dividing lines between regimes
        ref1, = plot.plot( [1, 1], [0, 200],linewidth = 1, linestyle = "--", c = "k")
        ref2, = plot.plot( [3.25, 3.25],[0, 200], linewidth = 1, linestyle = "--", c = "k")

    #if scale_height_version:
    #    # Reference lines for each scale height of separation
    #    for i in range(10):
    #        plot.plot(x_ref, y_ref + i + 1, linewidth = 1 + (i % 2), c = 'k', alpha = 0.4)

    #y_h4 = marginally_stable(x_ref, h = 0.04)
    #y_h6 = marginally_stable(x_ref, h = 0.06)
    #y_h8 = marginally_stable(x_ref, h = 0.08)
    #ref_h4, = plot.plot(x_ref, y_h4, linewidth = linewidth, c = colors[0], label = r"$h = 0.04$")
    #ref_h6, = plot.plot(x_ref, y_h6, linewidth = linewidth, c = colors[1], label = r"$h = 0.06$")
    #ref_h8, = plot.plot(x_ref, y_h8, linewidth = linewidth, c = colors[2], label = r"$h = 0.08$")

    # Actual Values (Initial, Secondary, Interior, Too Interior)
    i_h4 = plot.scatter(initial_density_maxima_h4 - initial_critical_maxima_h4, initial_critical_maxima_h4, s = 100 + size_offsets[initial_accretion_numbers_h4 - 1], c = colors[0], zorder = 100, alpha = alpha)
    i_h6 = plot.scatter(initial_density_maxima_h6 - initial_critical_maxima_h6, initial_critical_maxima_h6, s = 100 + size_offsets[initial_accretion_numbers_h6 - 1], c = colors[1], zorder = 100, alpha = alpha)
    i_h8 = plot.scatter(initial_density_maxima_h8 - initial_critical_maxima_h8, initial_critical_maxima_h8, s = 100 + size_offsets[initial_accretion_numbers_h8 - 1], c = colors[2], zorder = 100, alpha = alpha)

    s_h4 = plot.scatter(second_density_maxima_h4 - second_critical_maxima_h4, second_critical_maxima_h4, s = 175 + size_offsets[second_accretion_numbers_h4 - 1], c = colors[0], zorder = 100, marker = "*", alpha = alpha)
    s_h6 = plot.scatter(second_density_maxima_h6 - second_critical_maxima_h6, second_critical_maxima_h6, s = 175 + size_offsets[second_accretion_numbers_h6 - 1], c = colors[1], zorder = 100, marker = "*", alpha = alpha)
    #s_h8, = plot.scatter(second_critical_maxima_h8, second_critical_maxima_h8, s = 150, c = colors[2], zorder = 100, marker = "*")

    d_h4 = plot.scatter(interior_density_maxima_h4 - interior_critical_maxima_h4, interior_critical_maxima_h4, s = 100 + size_offsets[interior_accretion_numbers_h4 - 1], c = colors[0], zorder = 100, marker = "D", alpha = alpha)
    d_h6 = plot.scatter(interior_density_maxima_h6 - interior_critical_maxima_h6, interior_critical_maxima_h6, s = 100 + size_offsets[interior_accretion_numbers_h6 - 1], c = colors[1], zorder = 100, marker = "D", alpha = alpha)
    #d_h8, = plot.scatter(interior_critical_maxima_h8, interior_density_maxima_h8, s = 100, c = colors[2], zorder = 100, marker = "D")

    t_h4 = plot.scatter(too_interior_density_maxima_h4 - too_interior_critical_maxima_h4, too_interior_critical_maxima_h4, s = 100 + size_offsets[too_interior_accretion_numbers_h4 - 1], c = colors[0], zorder = 100, marker = "x", alpha = alpha)
    t_h6 = plot.scatter(too_interior_density_maxima_h6 - too_interior_critical_maxima_h6, too_interior_critical_maxima_h6, s = 100 + size_offsets[too_interior_accretion_numbers_h6 - 1], c = colors[1], zorder = 100, marker = "x", alpha = alpha)
    #d_h8, = plot.scatter(interior_critical_maxima_h8, interior_density_maxima_h8, s = 100, c = colors[2], zorder = 100, marker = "x")


    # Axes
    plot.xlim(1, 1.6)
    plot.ylim(0, 1.6)

    if scale_height_version:
        plot.ylim(0, 8)
        plot.xlim(0, 5)

    #plot.xscale("log"); plot.yscale("log")

    #plot.xticks([0.04, 0.06, 0.08, 0.1], ["0.04", "0.06", "0.08", "0.10"])
    #plot.yticks([0.1, 0.3, 1], ["0.1", "0.3", "1.0"])

    # Annotate
    if scale_height_version:
        unit1 = "(r_\mathrm{crit} - r_\mathrm{p}) / h"
        unit2 = "(r_\mathrm{pressure} - r_\mathrm{crit}) / h"
        plot.ylabel(r"$%s$" % unit1, fontsize = fontsize)
        #plot.ylabel(r"Density Maximum [$%s$]" % unit, fontsize = fontsize)
        plot.xlabel(r"Separation: $%s$" % unit2, fontsize = fontsize)
    else:
        unit = "r_\mathrm{p}"
        plot.xlabel(r"Radius of Critical Function Maximum [$%s$]" % unit, fontsize = fontsize)
        plot.ylabel(r"Radius of Density Maximum [%s$]" % unit, fontsize = fontsize)
    plot.title(r"Vortex Outcomes", y = 1.02, fontsize = fontsize + 2)

    #plot.text(0.0295, 2.35, "Unstable", fontsize = fontsize)
    #plot.text(0.078, 0.038, "Stable", fontsize = fontsize)
    
    #first_legend = plot.legend([t_h6, d_h6, s_h6, i_h6], ["No Vortex", "Interior Re-trigger", "Re-triggered Vortex", "Initial Vortex"], loc = "lower right", fontsize = fontsize - 2, scatterpoints = 1)
    first_legend = plot.legend([i_h6, s_h6, d_h6, t_h6], ["Initial Vortex", "Re-triggered Vortex", "Re-trigger (interior)", "No Vortex"], loc = "lower right", fontsize = fontsize - 2, scatterpoints = 1)
    second_legend = plot.legend([i_h4, i_h6, i_h8], [r"$h = 0.04$", r"$h = 0.06$", r"$h = 0.08$"], loc = "upper right", fontsize = fontsize - 2, scatterpoints = 1)
    ax = plot.gca().add_artist(first_legend)


    #if version is None:
    #    save_fn = "%s_onoGrid.png" % ()
    #else:
    #    save_fn = "v%04d_%s_onoGrid.png" % (version)
    if scale_height_version:
        plot.savefig("maximumGrid-h.png", bbox_inches = 'tight', dpi = dpi)
        plot.savefig("maximumGrid-h.pdf", bbox_inches = 'tight', dpi = dpi, format = "pdf")
    else:
        plot.savefig("maximumGrid.png", bbox_inches = 'tight', dpi = dpi)
        plot.savefig("maximumGrid.pdf", bbox_inches = 'tight', dpi = dpi, format = "pdf")

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = True)