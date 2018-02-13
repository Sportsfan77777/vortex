"""
plot 2-D density maps

python plotDensityMaps.py
python plotDensityMaps.py frame_number
python plotDensityMaps.py -1 <<<===== Plots a sample
python plotDensityMaps.py -m
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


# Beam Sizes for Extent Comparison
beam_sizes = np.array([1, 5, 10, 15, 20, 25, 30, 40])

def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = 2,
                         help = 'select two frames to compare intensity (first is for T = 10, second is for T = 1000)')

    # Directory Selection
    parser.add_argument('--dir1', dest = "directory1", default = '../taper10',
                         help = 'select first directory to compare intensity (first is for T = 10, second is for T = 10) (default: ../taper10)')
    parser.add_argument('--dir2', dest = "directory2", default = '../taper1000',
                         help = 'select second directory to compare intensity (first is for T = 10, second is for T = 1000) (default: ../taper1000)')

    parser.add_argument('-w', dest = "wavelength", type = float, default = 870,
                         help = 'wavelength (in um) (default: 870)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalIntensityComparison",
                         help = 'save directory (default: azimuthalIntensityComparison)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", nargs = '+', type = float, default = None,
                         help = 'max_y for each frame, or same for all (default: None)')
    parser.add_argument('-s', dest = "sliver_width", type = float, default = 2.0,
                         help = 'number of scale heights in sliver (default: 2.0)')

    parser.add_argument('-t', dest = "threshold", type = float, default = 0.5,
                         help = 'threshold for measuring extent (default: 0.5)')
    
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
default_directory = "taper1000/synthetic/lambda%04d/beam%03d" % (args.wavelength, args.beam_size)

fn = "../%s/id%04d_par.p" % (default_directory, args.id_number)
fargo_par = pickle.load(open(fn, "rb"))

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

wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

### Get Input Parameters ###

# Frames
frame_range = args.frames

# Directories
directory1 = "%s/synthetic" % args.directory1
directory2 = "%s/synthetic" % args.directory2

directories1 = ["%s/lambda%04d/beam%03d" % (directory1, wavelength, beam_size) for beam_size in beam_sizes]
directories2 = ["%s/lambda%04d/beam%03d" % (directory2, wavelength, beam_size) for beam_size in beam_sizes]

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
normalize = args.normalize

show = args.show
max_y = args.max_y

sliver_width = args.sliver_width

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

threshold = args.threshold

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

### Data ###

def get_extents(directories, frame):
    extents = np.zeros(len(directories))

    cwd = os.getcwd()

    for i, directory_i in enumerate(directories):
        # cd directory_i
        os.chdir(directory_i)

        # Get Data
        intensity_polar = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number, directory = "lambda%04d/beam%03d" % (args.wavelength, args.beam_size))
        extent = get_extent(intensity_polar, fargo_par, threshold = threshold, sliver_width = sliver_width):

        extents[i] = extent

        # cd back
        os.chdir(cwd)

    return extents


###############################################################################

##### PLOTTING #####

colors = ['#f20202', '#0609ef']

def make_plot(show = False):
    pass


##### Make Plots! #####

make_plot(show = show)

