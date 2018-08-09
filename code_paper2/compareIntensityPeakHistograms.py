"""
compare intensityPeakHistograms at different beam sizes
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
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import gridspec

from pylab import rcParams
from pylab import fromfile

import util
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')
    parser.add_argument('-b', dest = "beams", nargs = '+', type = int, default = [10, 20, 30, 40],
                         help = 'thresholds for marking edges of vortex with intensity (default: 0.5)')
    parser.add_argument('-t', dest = "thresholds", nargs = '+', type = float, default = [0.4, 0.5, 0.6, 0.7],
                         help = 'thresholds for marking edges of vortex with intensity (default: 0.5)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "intensityPeakHistograms",
                         help = 'save directory (default: intensityPeakHistograms)')

    # Save Data
    parser.add_argument('--save', dest = "save_data", action = 'store_true', default = False,
                         help = 'save data or not (default: do not save)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'max_y for plot (default: None)')

    parser.add_argument('--cum', dest = "cumulative", action = 'store_true', default = False,
                         help = 'normal or cumulative (default: not cumulative)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'labelsize of plot annotations (default: 16)')
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

### Get ID%04d Parameters ###
fn = "id%04d_par.p" % (args.id_number)
fargo_par = pickle.load(open("../beam010/%s" % fn, "rb"))

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
dust_surface_density_zero = surface_density_zero / 100
disk_mass = 2 * np.pi * dust_surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

planet_radius = fargo_par["Radius"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

beams = args.beams
thresholds = args.thresholds

# Save Data
save_data = args.save_data

# Plot Parameters (variable)
show = args.show
max_y = args.max_y
#if max_y is None:
#    max_y = 1

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

cumulative = args.cumulative

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

min_x = -90; max_x = 90
min_y = 0; max_y = 1

data = np.zeros((len(beams), len(frame_range)))
colors = ["b", "g", "y", "k"]

def make_plot(show = False):
    fig = plot.figure(figsize = (7, 5), dpi = dpi)
    ax = fig.add_subplot(111)

    # Get Data
    for i, (beam_i, threshold_i) in enumerate(zip(beams, thresholds)):
        data[i] = pickle.load(open("../beam%03d/id%04d_b%02d_t%02d_intensityPeaks.p" % (beam_i, id_number, beam_i, int(round(100.0 * threshold_i, 0))), "rb"))

    # Plot
    for i, beam_i in enumerate(beams):
        data_i = data[i]
        if cumulative:
            bins = np.linspace(min_x - 10, max_x + 10, 201) # Make this parameters
            hist = plot.hist(data_i, bins = bins, normed = True, color = colors[i], histtype = 'step', linewidth = linewidth, label = "%d" % beam_i, cumulative = True)
        else:
            bins = np.linspace(min_x - 10, max_x + 10, 21) # Make this parameters
            hist = plot.hist(data_i, bins = bins, normed = True, color = colors[i], histtype = 'step', linewidth = linewidth, label = "%d" % beam_i)

    # Minor Guidelines
    vertical = np.linspace(-30, 30, 7)
    horizontal = np.linspace(0, 1, 11)

    for vertical_i in vertical:
        plot.plot([vertical_i, vertical_i], [min_y, max_y], c = "k", linestyle = "--", alpha = alpha)
    #for horizontal_i in horizontal:
    #    plot.plot([min_x, max_x], [horizontal_i, horizontal_i], c = "k", alpha = alpha)

    # Axes
    plot.xlim(min_x, max_x)
    plot.ylim(min_y, max_y)

    xticks = np.linspace(min_x, max_x, 7)
    plot.xticks(xticks)

    plot.xlabel("Peak Offsets", fontsize = fontsize)
    plot.ylabel("Frequency", fontsize = fontsize)
    #plot.title("")

    # Legend
    plot.legend(loc = "upper left")

    # Save, Show, and Close
    frame_str = ""
    for i, frame_i in enumerate(args.frames):
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '-'

    if version is None:
        save_fn = "%s/intensityPeakHistograms_%s.png" % (save_directory, frame_str)
    else:
        save_fn = "%s/v%04d_intensityPeakHistograms_%s.png" % (save_directory, version, frame_str)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)



##### Make Plots! #####

make_plot(show = show)





