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
from readTitle import readTitle

from advanced import Parameters
from reader import Fields

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

master_directories = {}
#master_directories[67] = ["h06_nu7_a40_s5787-K60", "h06_nu7_a40_s3472-K60", "h06_nu7_a40_s2315-K60", "h06_nu7_a40_s1736-K60", "h06_nu7_a40_s1157-K60", "h06_nu7_a40_s0926-K60", "h06_nu7_a40_s0694-K60", "h06_nu7_a40_s0578-K60", "h06_nu7_a40_s0347-K60"]
master_directories[67] = ["h06_nu7_a40_s5787-K60", "h06_nu7_a40_s3472-K60", "h06_nu7_a40_s2315-K60", "h06_nu7_a40_s1736-K60", "h06_nu7_a40_s1157-K60", "h06_nu7_a40_s0694-K60", "h06_nu7_a40_s0347-K60"]

master_disc_masses = {}
#master_disc_masses[67] = [5.0, 3.0, 2.0, 1.5, 1.0, 0.8, 0.6, 0.5, 0.3]
master_disc_masses[67] = [5.0, 3.0, 2.0, 1.5, 1.0, 0.6, 0.3]

master_start_times = {}
#master_start_times[67] = [10, 10, 310, 270, 430, 530, 750, 980, -1]
master_start_times[67] = [10, 10, 310, 270, 430, 750, -1]

master_end_times = {}
#master_end_times[67] = [-1, -1, -1, -1, 3350, 1700, 1920, 1880, -1]
master_end_times[67] = [-1, -1, -1, -1, 3350, 1920, -1]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('choice', type = int,
                         help = 'choice of directories')
    parser.add_argument('--dir', dest = "save_directory", default = "mass",
                         help = 'save directory (default: mass)')

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
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 22,
                         help = 'fontsize of plot annotations (default: 22)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'fontsize of plot annotations (default: 3)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

directories = master_directories[args.choice]
disc_masses = master_disc_masses[args.choice]
start_times = master_start_times[args.choice]
end_times = master_end_times[args.choice]

### Get Fargo Parameters ###
p = Parameters(directory = "../" + directories[0])

num_rad = p.ny; num_theta = p.nx
r_min = p.ymin; r_max = p.ymax

surface_density_zero = p.sigma0

taper_time = p.masstaper

viscosity = p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters(directory = "../" + directories[0])

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
    x_min = 0; x_max = 5000
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

##### PLOTTING #####

colors = ['k', 'darkviolet', 'cornflowerblue', 'sandybrown', 'darkorange', 'r', 'grey']
labelsize = 19
size = 100

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 12), dpi = dpi)
    plot.subplots_adjust(hspace = 0.12)

    ######### TOP PANEL #########
    ax = fig.add_subplot(211)

    # Iterate
    for i, directory in enumerate(directories):
        # Label
        scale_height = float(directories[0].split("_")[0][1:]) / 100.0
        log_viscosity = float(directories[0].split("_")[1][2:]) - 2.0
        disc_mass = disc_masses[i]

        start_time = start_times[i]
        end_time = end_times[i]

        #label = r"$h =$ $%.02f$, $\alpha_\mathrm{visc} = 3 \times 10^{-%d}$, A = %.02f" % (scale_height, log_viscosity, accretion_rate)
        #label = r"$A = %.02f$" % (accretion_rate)
        label = disc_mass

        legend_text = r"$\Sigma_\mathrm{0}$ $/$ $\Sigma_\mathrm{base}$"

        # Data
        data = np.loadtxt("../%s/planet0.dat" % directory)
        times = data[:, 0]
        base_mass = data[:, 7]
        accreted_mass = data[:, 8]

        total_mass = base_mass + accreted_mass

        if negative:
            negative_mass = data[:, 13]
            total_mass -= negative_mass

        ### Plot ###
        # Basic
        x = times
        y = total_mass / jupiter_mass
        result = plot.plot(x, y, c = colors[i % len(colors)], linewidth = linewidth - 1, zorder = 99, label = label)

        # Vortex Lifetime
        if start_time > 0:
            start_time_i = az.my_searchsorted(x, start_time)
            if end_time > 0:
                end_time_i = az.my_searchsorted(x, end_time)
            else:
                end_time_i = -1

            result = plot.plot(x[start_time_i:end_time_i], y[start_time_i:end_time_i], c = colors[i % len(colors)], linewidth = linewidth + 3, zorder = 99)

            plot.scatter(x[start_time_i], y[start_time_i], c = colors[i % len(colors)], s = 150, marker = "o", zorder = 120)

            if end_time > 0:
                plot.scatter(x[end_time_i], y[end_time_i], c = colors[i % len(colors)], s = 175, marker = "H", zorder = 120)

    legend = plot.legend(loc = "upper right", fontsize = fontsize - 4, title = legend_text, facecolor = 'white', framealpha = 0.9)
    legend.set_zorder(150)
    ax.get_legend().get_title().set_fontsize(fontsize - 4)

    # Axes
    if args.max_y is None:
        x_min_i = np.searchsorted(x, x_min)
        x_max_i = np.searchsorted(x, x_max)
        max_y = 1.1 * max(y[x_min_i : x_max_i])
    else:
        max_y = args.max_y

    plot.xlim(x_min, x_max)
    plot.ylim(0, max_y)

    #title = readTitle()

    unit = "planet orbits"
    #plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel(r"$M_\mathrm{p}$ [$M_\mathrm{Jup}$]", fontsize = fontsize)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    title = r"$h = %.02f$          $\alpha \approx %s \times 10^{%d}$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2)
    #title = r"$h = %.02f$          $\alpha_\mathrm{disk} = 3 \times 10^{-%d}$" % (scale_height, log_viscosity)
    plot.title("%s" % (title), y = 1.015, fontsize = fontsize + 2)
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    ######### BOTTOM PANEL #########
    ax = fig.add_subplot(212)

    # Iterate
    for i, directory in enumerate(directories):
        # Label
        scale_height = float(directories[0].split("_")[0][1:]) / 100.0
        log_viscosity = float(directories[0].split("_")[1][2:]) - 2.0
        disc_mass = disc_masses[i]

        start_time = start_times[i]
        end_time = end_times[i]

        #label = r"$h =$ $%.02f$, $\alpha_\mathrm{visc} = 3 \times 10^{-%d}$, A = %.02f" % (scale_height, log_viscosity, accretion_rate)
        #label = r"$A = %.02f$" % (accretion_rate)
        label = disc_mass

        legend_text = r"$\Sigma_\mathrm{0}$ $/$ $\Sigma_\mathrm{base}$"

        # Data
        data = np.loadtxt("../%s/planet0.dat" % directory)
        times = data[:, 0]
        base_mass = data[:, 7]
        accreted_mass = data[:, 8]

        total_mass = base_mass + accreted_mass

        if negative:
            negative_mass = data[:, 13]
            total_mass -= negative_mass

        accretion = total_mass[1:] - total_mass[:-1]

        # Filter out repeats
        times = (times[1:])[accretion > 0]
        accretion = accretion[accretion > 0]

        ### Plot ###
        # Basic
        x = times[9:]
        y = accretion[9:] / jupiter_mass
        result = plot.plot(x, y, c = colors[i % len(colors)], linewidth = linewidth - 2, zorder = 99, label = label)

        # Vortex Lifetime
        if start_time > 0:
            start_time_i = az.my_searchsorted(x, start_time)

            if end_time > 0:
                end_time_i = az.my_searchsorted(x, end_time)
            else:
                end_time_i = -1

            result = plot.plot(x[start_time_i:end_time_i], y[start_time_i:end_time_i], c = colors[i % len(colors)], linewidth = linewidth + 1, zorder = 99)

            plot.scatter(x[start_time_i], y[start_time_i], c = colors[i % len(colors)], s = 150, marker = "o", zorder = 120)
            if end_time > 0:
                plot.scatter(x[end_time_i], y[end_time_i], c = colors[i % len(colors)], s = 175, marker = "H", zorder = 120)

    #plot.legend(loc = "upper right", fontsize = fontsize - 4)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(10**(-6), 10**(-2))

    plot.yscale("log")

    #title = readTitle()

    unit = "planet orbits"
    plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel(r"$\dot{M}_\mathrm{p}$ [$M_\mathrm{Jup}$ $/$ $T_\mathrm{p}$]", fontsize = fontsize)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/severalMassAndAccretionRates_choice%d.png" % (save_directory, args.choice)
    else:
        save_fn = "%s/v%04d_severalMassAndAccretionRates_choice%d.png" % (save_directory, version, arg.choice)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
