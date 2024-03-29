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

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "averagedDensity",
                         help = 'save directory (default: gasDensityMaps)')
    parser.add_argument('--mpi', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpio output files (default: use dat)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'maximum density (default: 1.1 times the max)')

    parser.add_argument('--zero', dest = "zero", action = 'store_true', default = False,
                         help = 'plot density at t = 0 for reference (default: do not do it!)')

    parser.add_argument('--data', dest = "data", default = None,
                         help = 'compare to data from another directory (default: do not do it!)')
    parser.add_argument('--data2', dest = "data2", default = None,
                         help = 'compare to data from another directory two (default: do not do it!)')
    parser.add_argument('--data3', dest = "data3", default = None,
                         help = 'compare to data from another directory three (default: do not do it!)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
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

viscosity = p.nu
scale_height = p.aspectratio

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

"""
fargo_par = util.get_pickled_parameters()

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

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

mpi = args.mpi

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

### Add new parameters to dictionary ###
#fargo_par["rad"] = rad
#fargo_par["theta"] = theta

###############################################################################

##### PLOTTING #####

def make_plot(show = False):
    # Figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)

    # Curves
    frame_range = pickle.load(open("./excess_mass_frames.p", "rb"))
    mass_over_time = pickle.load(open("./excess_mass_values.p", "rb"))
    plot.plot(frame_range, mass_over_time, linewidth = linewidth)

    if args.data:
        frame_range_data = pickle.load(open("%s/excess_mass_frames.p" % args.data, "rb"))
        mass_over_time_data = pickle.load(open("%s/excess_mass_values.p" % args.data, "rb"))
        plot.plot(frame_range_data, mass_over_time_data, linewidth = linewidth, label = "data")

    if args.data2:
        frame_range_data2 = pickle.load(open("%s/excess_mass_frames.p" % args.data2, "rb"))
        mass_over_time_data2 = pickle.load(open("%s/excess_mass_values.p" % args.data2, "rb"))
        plot.plot(frame_range_data2, mass_over_time_data2, linewidth = linewidth, label = "data2")

    if args.data3:
        frame_range_data3 = pickle.load(open("%s/excess_mass_frames.p" % args.data3, "rb"))
        mass_over_time_data3 = pickle.load(open("%s/excess_mass_values.p" % args.data3, "rb"))
        plot.plot(frame_range_data3, mass_over_time_data3, linewidth = linewidth, label = "data3")

    # Annotate
    #this_title = readTitle()
    title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Excess Mass", fontsize = fontsize)
    plot.title(title1, fontsize = fontsize)

    #plot.legend(loc = "upper left")

    # Limits
    plot.xlim(0, frame_range[-1])
    #plot.ylim(0.0, 1.0)

    plot.yscale('log')

    # Save + Close
    plot.savefig("excessMass.png")
    plot.show()

    plot.close()

##### Make Plots! #####

make_plot(show = show)
