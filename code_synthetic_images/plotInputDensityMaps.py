"""
plot 2-D density maps

python plotInputDensityMaps.py
python plotInputDensityMaps.py frame_number
python plotInputDensityMaps.py -1 <<<===== Plots a sample
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

import matplotlib
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util

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
    parser.add_argument('-s', dest = "run_parallel", action = 'store_false', default = True,
                         help = 'for multiple plots, run in serial (default: run in parallel)') 
    ######## ^^^^^^^ Switch run_parallel with number of cores!!! ##########

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "gasDensityMaps",
                         help = 'save directory (default: gasDensityMaps)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')

    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of plot parameters (default: None)')
    parser.add_argument('--r_range', dest = "r_lim", type = int, nargs = 2, default = None,
                         help = 'id number for this set of plot parameters (default: [r_min, r_max])')
    

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "viridis",
                         help = 'color map (default: viridis)')
    parser.add_argument('--cmax', dest = "cmax", type = int, default = 2.5,
                         help = 'maximum density in colorbar (default: 2.5)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get FARGO Parameters ###
fargo_par = util.get_pickled_parameters(directory = "../cm-size") # Retrieve parameters from *.par file

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
surface_density_zero = fargo_par["Sigma0"]
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
alpha_viscosity = fargo_par["AlphaViscosity"]

### Get Input Parameters ###

# Frames
if len(args.frames) == 1:
    frame_range = args.frames
elif len(args.frames) == 3:
    start = args.frames[0]; end = args.frames[1]; rate = args.frames[2]
    frame_range = range(start, end + 1, rate)

    run_parallel = args.run_parallel
else:
    print "Error: Must supply 1 or 3 frame arguments\nWith one argument, plots single frame\nWith three arguments, plots range(start, end + 1, rate)"
    exit()

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show

id_number = args.id_number
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]

# Plot Parameters (constant)
cmap = args.cmap
clim = [0, args.cmax]

fontsize = args.fontsize
dpi = args.dpi

###############################################################################

##### PLOTTING #####

def make_plot(frame, show = False):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')

    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    density = (fromfile("i%04d_gasddens%d.dat" % (id_number, frame)).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    # Center Data
    middle = np.searchsorted(theta, np.pi)
    normalized_density = np.roll(normalized_density, middle, axis = -1)

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi)
    result = ax.pcolormesh(x, y, np.transpose(normalized_density), cmap = cmap)

    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)
    plot.title("Gas Density Map (t = %.1f)" % (orbit), fontsize = fontsize + 1)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)

    # Save, Show,  and Close
    if id_number is None:
        save_fn = "%s/densityMap_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/id%04d_densityMap_%04d.png" % (save_directory, id_number, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

if len(frame_range) == 1:
    make_plot(frame_range[0], show = show)
else:
    if run_parallel:
        p = Pool(4) # default number of processes is multiprocessing.cpu_count()
        p.map(make_plot, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            make_plot(frame)

