"""
tracks total dust mass in annulus encapsulating the vortex
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams

import util
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot total dust mass in vortex annulus over time."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select range(start, end, rate)')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = ".",
                         help = 'save directory (default: current directory)')
    parser.add_argument('--name', dest = "name", default = "default",
                         help = 'save name suffix (default: default)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')

    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'radial range in plot (default: None)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 18,
                         help = 'fontsize of plot annotations (default: 18)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 15,
                         help = 'labelsize of plot annotations (default: 15)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'linewidths in plot (default: 3)')
    parser.add_argument('--alpha', dest = "alpha", type = float, default = 0.65,
                         help = 'line transparency in plot (default: 0.65)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
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

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist
name = args.name

# Plot Parameters (variable)
show = args.show
max_y = args.max_y

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Helper Functions ###
def get_total_mass(frame, num_scale_heights = 3):
    # Get Data
    density = util.read_dust_data(frame, fargo_par) # Dust!!!!

    # Find Center
    avgDensity = np.average(density, axis = -1)
    peak_rad, peak_density = az.get_radial_peak(avgDensity, fargo_par)

    # Zoom in on annulus
    start_rad = peak_rad - 0.5 * num_scale_heights * scale_height
    end_rad = peak_rad + 0.5 * num_scale_heights * scale_height

    start_rad_i = np.searchsorted(rad, start_rad)
    end_rad_i = np.searchsorted(rad, end_rad)

    rad_annulus = rad[start_rad_i : end_rad_i]
    density_annulus = density[start_rad_i : end_rad_i]

    # Multiply by Grid Cell size
    dr = rad[1] - rad[0]; dphi = theta[1] - theta[0]
    density_annulus = (dr * dphi) * rad_annulus[:, None] * density_annulus # (r * dr * d\phi) * \rho

    total_mass = np.sum(density_annulus)

    return total_mass

###############################################################################

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Line Plot ###
    x = np.array(frame_range)
    y = np.array([get_total_mass(frame) for frame in x])
    plot.plot(x, y, linewidth = linewidth)

    # Annotate Axes
    plot.xlabel("Time (planet orbits)", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)
    plot.title("Total Dust Mass in Vortex Annulus")

    #plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 1.0)) # outside of plot

    # Axes
    plot.xlim(frame_range[0], frame_range[-1])
    if max_y is not None:
        plot.ylim(0, max_y)

    # Save, Show, and Close
    save_fn = "totalVortexDustMass_%s.png" % name
    plot.savefig(save_fn, bbox_inches = 'tight')

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)