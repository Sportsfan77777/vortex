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
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
import azimuthal as az
#from readTitle import readTitle

from advanced import Parameters
#from reader import Fields

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "radii",
                         help = 'save directory (default: mass)')

    # Reference
    parser.add_argument('--ref', dest = "ref", type = int, default = 0,
                         help = 'reference taper time for prescribed growth curve (default: no reference)')
    parser.add_argument('--compare', dest = "compare", nargs = '+', default = ["../h10_nu5_b0-Mt3-768x768-case4-stockholm-migrating", "../h10_nu5_b45-Mt3-768x768-case4-stockholm-migrating-no-extra-gap-torque", "../h10_nu5_b45-Mt3-768x768-case12-stockholm-migrating"],
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

    parser.add_argument('--offset', dest = "offset", type = float, default = 0.0,
                         help = 'time offset for compare (default: 0.0)')

    parser.add_argument('--negative', dest = "negative", action = 'store_true', default = False,
                         help = 'add negative mass (default: do not)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 20,
                         help = 'fontsize of plot annotations (default: 20)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'fontsize of plot annotations (default: 3)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
p = Parameters()

num_rad = p.ny; num_theta = p.nx
r_min = p.ymin; r_max = p.ymax

surface_density_zero = p.sigma0

taper_time = p.masstaper

scale_height = p.aspectratio
viscosity = p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]


"""
num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]


taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

size = fargo_par["PSIZE"]
"""

### Get Input Parameters ###

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Reference
ref = args.ref

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
if args.r_lim is None:
    x_min = 0; x_max = 1000
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
max_y = args.max_y

negative = args.negative

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

### Add new parameters to dictionary ###
#fargo_par["rad"] = rad
#fargo_par["theta"] = theta

###############################################################################

labelsize = 18
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    data = np.loadtxt("planet0.dat")
    times = data[:, 0]
    planet_x = data[:, 1]
    planet_y = data[:, 2]
    planet_radii = np.sqrt(np.power(planet_x, 2) + np.power(planet_y, 2))

    cwd = os.getcwd().split("/")[-1]

    ### Plot ###
    x = times
    y = planet_radii
    result = plot.plot(x, y, linewidth = linewidth + 1, zorder = 199, label = "Full Model")

    if args.ref > 0:
        x = times
        y_ref = np.power(np.sin((np.pi / 2) * (1.0 * times / args.ref)), 2) * 1.0
        plot.plot(x, y_ref, linewidth = linewidth, alpha = 0.5)

    #colors = ["k", "C1", "gold", "indigo", "C4", "grey"]
    colors = ["k", "gold", "indigo"]
    #labels = ["Just Viscosity", "FULL (No Viscosity)",  "FULL (No Extra Gap Torque)", "K20 (Viscous)", "K20 (No Viscosity)", "Just Inviscid"]
    labels = ["Just Viscosity",  "Full (No Extra Gap Torque)", "K20"]
    dashes = ["-.", "-", "-"]
    #orders = [50, 0, 0, 0, 0, 50]
    orders = [50, 0, 0]
    #alphas = [0.5, 1, 1, 1, 1, 0.5]
    alphas = [0.5, 1, 1]
    if args.compare is not None:
        for i, directory in enumerate(args.compare):
            data_comp = np.loadtxt("%s/planet0.dat" % directory)
            times = data_comp[:, 0] + args.offset

            planet_x = data_comp[:, 1]
            planet_y = data_comp[:, 2]
            planet_radii = np.sqrt(np.power(planet_x, 2) + np.power(planet_y, 2))

            x_comp = times
            y_comp = planet_radii
            result = plot.plot(x_comp, y_comp, c = colors[i], linewidth = linewidth, zorder = 99 - i - orders[i], label = labels[i], linestyle = dashes[i], alpha = alphas[i])

        h, l = ax.get_legend_handles_labels()
        handles_r = h[::-1]; labels_r = l[::-1]
        plot.legend(loc = "lower left", fontsize = fontsize - 5, handles = handles_r, labels = labels_r)

    # Axes
    if args.max_y is None:
        x_min_i = np.searchsorted(x, x_min)
        x_max_i = np.searchsorted(x, x_max)
        max_y = 1.1 * max(y[x_min_i : x_max_i])
    else:
        max_y = args.max_y

    plot.xlim(x_min, x_max)
    plot.ylim(0.3, max_y)

    #title = readTitle()

    unit = "orbits"
    plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel(r"Planet semi-major axis [$r_p$]", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    #title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)
    title1 = r"$h = %.2f$    $M_\mathrm{p} = %.2f\ M_\mathrm{Jup}$" % (scale_height, planet_mass)
    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title1), y = 1.015, fontsize = fontsize + 1)
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    text_nu = r"$\nu = 10^{%d}$" % (int(np.log(viscosity) / np.log(10)))
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    plot.text(0.95 * x_range + x_min, 0.9 * (plot.ylim()[-1] - plot.ylim()[0]) + plot.ylim()[0], text_nu, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')


    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1]

    if version is None:
        save_fn = "%s/%s_radiiOverTime.png" % (save_directory, directory_name)
    else:
        save_fn = "%s/v%04d_%s_radiiOverTime.png" % (save_directory, version, directory_name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
