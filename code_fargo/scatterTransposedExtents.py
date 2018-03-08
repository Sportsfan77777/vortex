"""
scatter plot of azimuthal extents

Usage:
python scatterExtents.py
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
from matplotlib import patches as patch

from itertools import groupby

from pylab import rcParams # replace with rc ???
from pylab import fromfile

# Dictionaries
cases = [1, 10, 25, 250]
tapers = {}
extents = {}

# M = 1, v = 10^-6
case = cases[0]
#tapers[case] = [10, 50, 100, 150, 250, 500, 1000, 3000]
#extents[case] = [180, 180, 180, 240, 240, 240, 400, 800]
tapers[case] = [10, 50, 100, 150, 250, 500, 1000, 3000]
extents[case] = [180, 180, 180, 240, 240, 240, 400, 800]

# M = 1, v = 10^-7
case = cases[1]
#tapers[case] = [10, 100, 150, 250, 500, 1000, 2000]
#extents[case] = [120, 120, 120, 240, 240, 240, 240]
tapers[case] = [10, 50, 100, 150, 250, 500, 1000, 2000]
extents[case] = [120, 120, 120, 120, 240, 240, 240, 240]

# M = 5, v = 10^-6
case = cases[2]
#tapers[case] = [2, 100, 200, 400]
#extents[case] = [90, 120, 180, 240]
tapers[case] = [2, 50, 100, 150, 200, 400]
extents[case] = [90, 120, 120, 120, 180, 240]

# M = 5, v = 10^-7
case = cases[3]
#tapers[case] = [2, 200, 400, 800]
#extents[case] = [120, 120, 180, 180]
tapers[case] = [2, 50, 100, 150, 200, 400, 800]
extents[case] = [120, 120, 120, 120, 120, 180, 180]


## Helper Functions ##
x_min = np.log(30); x_max = np.log(5000)
y_min = np.log(0.5); y_max = np.log(600)

def map_to_y(y_cor):
    offset = (x_max - x_min) / (y_max - y_min)
    intercept = (x_min - y_min)
    return offset * y_cor + intercept

### PLOTTING ###
fontsize = 17
labelsize = 14
dpi = 100

colors = ["gold", "k", "b"]

e_base = 0.25
e_width_c = 1.0 * e_base # concentrated
e_width_i = 1.0 * e_base # intermediate
e_width_e = 1.0 * e_base # elongated
e_height_c = 1.0 * e_base # concentrated
e_height_i = 2.0 * e_base # intermediate
e_height_e = 4.0 * e_base # elongated

def make_plot(show = True):
    # Make Figure
    fig = plot.figure(figsize = (6, 6), dpi = dpi)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')

    # Scatter Points
    for c, case in enumerate(cases):
        tapers_i = tapers[case]
        extents_i = extents[case]

        for i, (taper, extent) in enumerate(zip(tapers_i, extents_i)):
            cutoff_i = int(sys.argv[1]) # plot up to cutoff for each case

            if extent > 360:
                # No Vortex!
                if i <= cutoff_i:
                    plot.scatter([np.log(taper)], [map_to_y(np.log(case))], marker = "x", c = "r", s = 400)
            else:
                # Vortex!
                if extent > 180:
                    # Elongated!
                    extent_type = "x"
                    e_width = e_width_e; e_height = e_height_e
                    color = colors[2]
                    # For SVM
                    #elongated_x.append(np.log(case) / np.log(10))
                    #elongated_y.append(np.log(taper) / np.log(10))
                elif extent == 180:
                    # Intermediate!
                    extent_type = "_"
                    e_width = e_width_i; e_height = e_height_i
                    color = colors[1]
                else:
                    # Concentrated!
                    extent_type = "o"
                    e_width = e_width_c; e_height = e_height_c
                    color = colors[0]
                    # For SVM
                    #concentrated_x.append(np.log(case) / np.log(10))
                    #concentrated_y.append(np.log(taper) / np.log(10))

                if i <= cutoff_i:
                    ellipse = patch.Ellipse(xy = (np.log(taper), map_to_y(np.log(case))), width = e_width, height = e_height, fc = color, zorder = 100)
                    fig.gca().add_artist(ellipse)

    # Extra Lines
    split_y = map_to_y(np.log(16))
    plot.plot([0.62 * x_min, x_max], [split_y, split_y], c = 'k', clip_on = False)

    # Axes
    #plot.xlim(30, 3000); plot.ylim(0.5, 400)
    #plot.xlim(np.log(30), np.log(3000)); plot.ylim(np.log(0.5), np.log(400))
    ax.set_xlim(x_min, x_max); ax.set_ylim(x_min, x_max)

    #plot.xscale("log"); plot.yscale("log")

    x_axis_labels = ["50", "100", "200", "500", "1000", "2000"]
    x_tick_locations = np.array([np.log(int(x)) for x in x_axis_labels])

    ax.set_xticks(x_tick_locations)
    ax.set_xticklabels(x_axis_labels, fontsize = labelsize)

    y_axis_labels = ["M1v4", "M1v5", "M5v4", "M5v5"]
    y_axis_labels = ["-4", "-5", "-4", "-5"]

    y_base = r"   $3 \times 10^{%d}$"
    y_base = r"$\alpha = 3 \times 10^{%d}$"
    y_axis_labels = [y_base % -4, y_base % -5, y_base % -4, y_base % -5]
    y_tick_locations = [map_to_y(np.log(c)) for c in cases]

    ax.set_yticks(y_tick_locations)
    ax.set_yticklabels(y_axis_labels, fontsize = labelsize + 1)

    ax2 = ax.twiny()
    x2_axis_labels = ["250", "500", "1000", "2500", "5000", ""]
    x2_tick_locations = np.array([np.log(int(x) * 0.94) for x in x_axis_labels])

    ax2.set_xlim(x_min, x_max)
    ax2.set_xticks(x2_tick_locations)
    ax2.set_xticklabels(x2_axis_labels, fontsize = labelsize)

    # Annotate Axes
    ax.set_xlabel(r"Time to grow to $1\ M_\mathrm{Jup}$ (planet orbits)", fontsize = fontsize)
    ax2.set_xlabel(r"Time to grow to $5\ M_\mathrm{Jup}$ (planet orbits)", fontsize = fontsize, labelpad = 12)
    #plot.ylabel(r"$M_\mathrm{p}^2 / \alpha_\mathrm{disk}$", fontsize = fontsize)
    #plot.title("Azimuthal Extents", y = 1.01, fontsize = fontsize + 2)

    plot.text(0.75 * x_min, 1.00 * x_max, "Stronger\nVortices", horizontalalignment = 'center', fontsize = fontsize - 1)
    plot.text(0.75 * x_min, 0.95 * x_min, "Weaker\nVortices", horizontalalignment = 'center', fontsize = fontsize - 1)

    # Annotate Masses
    mass_x = np.log(40)
    plot.text(mass_x, map_to_y(np.log(75)), r"$5\ M_\mathrm{Jup}$", horizontalalignment = 'left', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 1)
    plot.text(mass_x, map_to_y(np.log(3)), r"$1\ M_\mathrm{Jup}$", horizontalalignment = 'left', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 1)

    # Save, Show, and Close
    save_fn = "v5_transposed_extents_%d.png" % cutoff_i
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        pass
        #plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)



make_plot()