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
    parser.add_argument('--dir', dest = "save_directory", default = "extentsByBeam",
                         help = 'save directory (default: extentsByBeam)')

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
default_beam_size = 10
default_directory = "taper1000/synthetic/lambda%04d/beam%03d" % (args.wavelength, default_beam_size)

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

    for i, directory_i in enumerate(directories):
        # Get Data
        intensity_polar = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number, directory = directory_i)
        extent = az.get_extent(intensity_polar, fargo_par, normalize = True, threshold = threshold, sliver_width = sliver_width)

        extents[i] = extent * (180.0 / np.pi)

    return extents


###############################################################################

##### PLOTTING #####

colors = ['#f20202', '#0609ef']
labels = [r"$T_\mathrm{growth} = 10$", r"$T_\mathrm{growth} = 1000$"]

def make_plot(show = False):
    fig = plot.figure(figsize = (7, 5), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    extents1 = get_extents(directories1, frame_range[0])
    extents2 = get_extents(directories2, frame_range[1])

    extent_arrays = [extents1, extents2]

    # Plot
    x = beam_sizes / planet_radius
    x_arc = 1.0 * beam_sizes / distance

    for i, extent_array in enumerate(extent_arrays):
        plot.plot(x_arc, extent_array, c = colors[i], linewidth = linewidth, label = labels[i])

    difference = extents2 - extents1
    plot.plot(x_arc, difference, c = "k", linewidth = linewidth - 1, linestyle = "--", label = r"$\mathrm{Difference}$")

    # Axes
    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)
    plot.ylim(0, 360)

    # Annotate Axes
    plot.xlabel(r"$\mathrm{Beam\ Diameter}$ [$^{\prime \prime}$]", fontsize = fontsize + 2)
    plot.ylabel(r"$\phi_\mathrm{extent}$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)

    plot.legend(loc = "upper left")

    # Title
    title = r"Azimuthal Extents at $I / I_0 = %.1f$" % (threshold)
    plot.title("%s" % (title), y = 1.10, fontsize = fontsize + 3)

    # Second x-axis
    xlim = plot.xlim()
    ax2 = ax.twiny()
    ax2.set_xlim(xlim[0] * (1.0 * distance / planet_radius), xlim[-1] * (1.0 * distance / planet_radius))

    # Save, Show, and Close
    frame_str = ""
    for i, frame_i in enumerate(frame_range):
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    png = "png"; pdf = "pdf"
    if version is None:
        save_fn = "%s/extentsByBeam_%s.%s" % (save_directory, frame_str, png)
        pdf_save_fn = "%s/extentsByBeam__%s.%s" % (save_directory, frame_str, pdf)
    else:
        save_fn = "%s/v%04d_extentsByBeam_%s.%s" % (save_directory, version, frame_str, png)
        pdf_save_fn = "%s/v%04d_extentsByBeam_%s.%s" % (save_directory, version, frame_str, pdf)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)
    plot.savefig(pdf_save_fn, bbox_inches = 'tight', format = "pdf")

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)

