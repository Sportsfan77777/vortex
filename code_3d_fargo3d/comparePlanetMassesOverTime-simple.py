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
import argparse

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

###############################################################################

file_prefixes = ["h06_nu7_a40_s1157-K60-2D", "h06_nu7_a40_s1157-K60-3D"]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = ".",
                         help = 'save directory (default: .)')

    # Reference
    parser.add_argument('--ref', dest = "ref", type = int, default = 0,
                         help = 'reference taper time for prescribed growth curve (default: no reference)')
    parser.add_argument('--compare', dest = "compare", nargs = '+', default = None,
                         help = 'select directories to compare planet growth rates')


    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'maximum density (default: 1.1 times the max)')

    parser.add_argument('--negative', dest = "negative", action = 'store_true', default = False,
                         help = 'add negative mass (default: do not)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 20,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 5,
                         help = 'fontsize of plot annotations (default: 5)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Reference
ref = args.ref

# Plot Parameters (variable)
show = args.show

version = args.version
if args.r_lim is None:
    x_min = 0; x_max = 3000
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
max_y = args.max_y

negative = args.negative

jupiter_mass = 1e-3

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

### Add new parameters to dictionary ###
#fargo_par["rad"] = rad
#fargo_par["theta"] = theta

###############################################################################

##### PLOTTING #####
labelsize = 17
alpha = 0.8

#colors = ["b", "gold", "#17becf", "orange"]
colors = ["b", "gold"]
#labels = [r"2-D $256 \times 256$", r"3-D $256 \times 256$", r"2-D $512 \times 512$", r"3-D $512 \times 512$"]
labels = [r"2-D  $256 \times 256$", r"3-D  $256 \times 256$"]

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    if negative:
        negative_mass = data[:, 13]
        total_mass -= negative_mass

    ### Plot ###
    max_y = 0.0
    for i, prefix in enumerate(file_prefixes):
        # Parameters
        surface_density_zero = float(prefix.split("_")[3][1:5]) * 1e-7
        scale_height = float(prefix.split("_")[0][1:]) / 100.0
        viscosity = float(prefix.split("_")[1][2:]) - 2.0

        alpha_coefficent = "3"
        if scale_height == 0.08:
            alpha_coefficent = "1.5"
        elif scale_height == 0.04:
            alpha_coefficent = "6"

        # Data
        data = np.loadtxt("%s-planet0.dat" % prefix)
        times = data[:, 0]
        base_mass = data[:, 7]
        accreted_mass = data[:, 8]
        total_mass = base_mass + accreted_mass

        x = times
        y = total_mass / jupiter_mass
        result = plot.plot(x, y, c = colors[i], linewidth = linewidth, zorder = 99, label = labels[i], alpha = alpha)

        # Axis
        x_min_i = np.searchsorted(x, x_min)
        x_max_i = np.searchsorted(x, x_max)
        max_y_temp = 1.1 * max(y[x_min_i : x_max_i])
        if max_y_temp > max_y:
            max_y = max_y_temp

    plot.legend(loc = "lower right", fontsize = fontsize - 4)

    # Axes
    if args.max_y is not None:
        max_y = args.max_y

    plot.xlim(x_min, x_max)
    plot.ylim(0, max_y)

    #title = readTitle()

    unit = "planet orbits"
    plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel(r"$M_p$ [$M_J$]", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    #x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    #y_text = 1.14

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    x_shift = 0.05
    y_text = 1.03; y_shift = 0.08

    text1 = r"$h = %.2f$" % (scale_height)
    plot.text(x_min - x_shift * x_range, y_text * plot.ylim()[-1], text1, horizontalalignment = 'left', fontsize = fontsize)
    text2 = r"$\alpha \approx %s \times 10^{-%d}$" % (alpha_coefficent, viscosity)
    plot.text(x_max + x_shift * x_range, (y_text) * plot.ylim()[-1], text2, horizontalalignment = 'right', fontsize = fontsize)

    #title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)

    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    #plot.title("%s" % (title1), y = 1.015, fontsize = fontsize + 1)
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    surface_density_base = 1.157e-4
    #title = r"$\Sigma_0$ $/$ $\Sigma_\mathrm{base} = %.1f$    $\Sigma_0 = %.3e$" % (surface_density_zero / surface_density_base, surface_density_zero)
    title = r"$\Sigma_0$ $/$ $\Sigma_\mathrm{base} = %.1f$" % (surface_density_zero / surface_density_base)
    plot.title("%s" % (title), y = 1.04, fontsize = fontsize + 2, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))

    # Text
    #text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    #text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')


    # Save, Show, and Close
    current_directory = os.getcwd().split("/")[-1]
    if version is None:
        save_fn = "%s/massOverTime-%s.png" % (save_directory, file_prefixes[-1])
    else:
        save_fn = "%s/v%04d_massComparisonOverTime-%s.png" % (save_directory, version, file_prefixes[-1])
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
