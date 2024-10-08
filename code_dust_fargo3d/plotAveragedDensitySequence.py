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

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "averagedDensitySequences",
                         help = 'save directory (default: averagedDensitySequences)')
    parser.add_argument('--mpi', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpiio output files (default: do not use mpi)')
    parser.add_argument('--merge', dest = "merge", type = int, default = 0,
                         help = 'number of cores needed to merge data outputs (default: 0)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'maximum density (default: 1.1 times the max)')

    parser.add_argument('--start', dest = "start", type = int, default = None,
                         help = 'vortex lifetime start (default: None)')
    parser.add_argument('--end', dest = "end", type = int, default = None,
                         help = 'vortex lifetime end (default: None)')
    
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

planet_mass = 1.0
taper_time = p.masstaper

scale_height = p.aspectratio
viscosity = p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

"""
num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

size = fargo_par["PSIZE"]
"""

### Get Input Parameters ###

# Frames
frame_range = args.frames

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
max_y = args.max_y

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

# Planet File
# Data
data = np.loadtxt("planet0.dat")
times = data[:, 0]; base_mass = data[:, 7]
accreted_mass = data[:, 8] / jupiter_mass

planet_x = data[:, 1]
planet_y = data[:, 2]

### Add new parameters to dictionary ###
#fargo_par["rad"] = rad
#fargo_par["theta"] = theta

###############################################################################

##### PLOTTING #####

linestyles = ["-", "-"]
#colors = ['k', 'b', 'cornflowerblue', '#17becf', '#8c564b', 'darkorange', 'r', 'gold']
#colors = ['k', 'b', '#17becf', '#8c564b', 'darkorange', 'darkviolet']
colors = ['k', 'b', '#17becf', 'darkviolet']
colors = ['k', 'firebrick', 'darkorange', '#8c564b']

labelsize = 18
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

alpha = 0.8

def make_plot(frame_range, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    for i, frame in enumerate(frame_range):
        density = fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)
        averagedDensity = np.average(density, axis = 1)
        normalized_density = averagedDensity / surface_density_zero

        ### Plot ###
        x = rad
        y = normalized_density
        result = plot.plot(x, y, c = colors[i % len(colors)], linewidth = linewidth + (i != 0) + (i == 2), linestyle = linestyles[i % 2], zorder = 99, label = r"$t$ $=$ $%d$ $T_\mathrm{p}$" % frame, alpha = alpha)

        this_x = planet_x[frame]
        this_y = planet_y[frame]
        planet_location = np.sqrt(np.power(this_x, 2) + np.power(this_y, 2))

        plot.scatter([planet_location], [1.05 * 10**(-3)], c = colors[i % len(colors)], s = 75, alpha = 0.8, clip_on = False)

    plot.legend(loc = "lower right", fontsize = fontsize - 5)

    # Axes
    if args.max_y is None:
        x_min_i = np.searchsorted(x, x_min)
        x_max_i = np.searchsorted(x, x_max)
        max_y = 1.1 * max(y[x_min_i : x_max_i])
    else:
        max_y = args.max_y

    plot.xlim(x_min, x_max)
    #plot.ylim(0, max_y)

    plot.ylim(10**(-3), 3.0)
    plot.yscale("log")

    # Annotate Axes
    orbit = (dt / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    current_mass += accreted_mass[frame]

    #title = readTitle()

    unit = "r_\mathrm{p}"
    plot.xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
    plot.ylabel(r"$<$ $\Sigma$ $/$ $\Sigma_0$ $>$", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    x_shift = 0.35; y_shift = 0.08
    x_text = 0.68; y_text = 0.84; 
    x_text_header = 0.5; y_text_header = 1.12

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    text1 = r"$h = %.2f$" % (scale_height)
    plot.text(x_text, y_text + y_shift, text1, horizontalalignment = 'left', fontsize = fontsize - 1, transform = ax.transAxes)
    text2 = r"$\alpha \approx %s \times 10^{%d}$" % (alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2)
    plot.text(x_text, y_text, text2, horizontalalignment = 'left', fontsize = fontsize - 1, transform = ax.transAxes)

    surface_density_base = 1.157e-4
    final_frame = 5000
    if final_frame > len(accreted_mass):
        final_frame = len(accreted_mass) - 1
    final_planet_mass = planet_mass + accreted_mass[final_frame]
    #title = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    title = r"$\Sigma_0$ $/$ $\Sigma_\mathrm{base} = %.1f$    ($M_\mathrm{p} = %.2f$ $M_\mathrm{Jup}$)" % (surface_density_zero / surface_density_base, final_planet_mass)
    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title1 = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    #title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title), y = 1.015, fontsize = fontsize + 1)
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    header = "Migrating"
    plot.text(x_text_header, y_text_header, header, horizontalalignment = 'center', fontsize = fontsize + 2, transform = ax.transAxes)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')

    if args.start is not None:
        text_start = r"$t_\mathrm{start}$ $=$ $%d$ $T_\mathrm{p}$" % args.start
        plot.text(0.55, 0.0025, text_start, fontsize = fontsize - 4, color = 'black', horizontalalignment = 'left')
    if args.end is not None:
        text_end = r"$t_\mathrm{end}$  $=$ $%d$ $T_\mathrm{p}$" % args.end
        plot.text(0.55, 0.0015, text_end, fontsize = fontsize - 4, color = 'black', horizontalalignment = 'left')

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1].split("-")[0]
    frame_string = ""
    for frame in frame_range:
        frame_string += ("%04d-" % frame)
    frame_string = frame_string[:-1] # get rid of trailing '-'

    if version is None:
        save_fn = "%s/averagedDensity_%s_%s.png" % (save_directory, directory_name, frame_string)
    else:
        save_fn = "%s/v%04d_averagedDensity_%s_%s.png" % (save_directory, version, directory_name, frame_string)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

# Iterate through frames

make_plot(frame_range, show = show)
