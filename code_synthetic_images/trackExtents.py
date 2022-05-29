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
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
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

# Grain Sizes for Extent Comparison
sizes = np.array([1, 0.3, 0.1, 0.03]) #, 0.01, 0.0001])
grain_sizes = np.array(["cm", "hcm", "mm", "hmm"]) #, "hum", "um"])

# Reverse
sizes = sizes[::-1]
grain_sizes = grain_sizes[::-1]

def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "extentsOverTime",
                         help = 'save directory (default: extentsOverTime)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", nargs = '+', type = float, default = None,
                         help = 'max_y for each frame, or same for all (default: None)')
    parser.add_argument('-s', dest = "sliver_width", type = float, default = 1.0,
                         help = 'number of scale heights in sliver (default: 1.0)')

    parser.add_argument('-t', dest = "threshold", type = float, default = 0.2,
                         help = 'threshold for measuring extent (default: 0.2)')

    parser.add_argument('--compare', dest = "compare", action = 'store_true', default = False,
                         help = 'compare the elongated vortex extents to the concentrated ones at the same threshold (default: do not compare)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'labelsize of plot annotations (default: 16)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'linewidths in plot (default: 3)')
    parser.add_argument('--alpha', dest = "alpha", type = float, default = 0.35,
                         help = 'line transparency in plot (default: 0.35)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
fn = "id%04d_par.p" % (args.id_number)
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
viscosity = fargo_par["Nu"]

#size = fargo_par["PSIZE"]

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

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

compare = args.compare

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

def get_extent(args):
    # Extract args
    i, frame = args

    # Get data and measure extent
    intensity = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number)
    extent, azimuthal_profile = az.get_extent(intensity, fargo_par, normalize = True, threshold = threshold, sliver_width = sliver_width)

    # Convert to degrees
    extent *= (180.0 / np.pi)

    # Store
    extents[i] = extent

###############################################################################

# Data
extents = mp_array("f", len(frame_range))
pool_args = [(i, frame) for i, frame in enumerate(frame_range)]    

p = Pool(num_cores)
p.map(get_extent, pool_args)
p.terminate()


##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

colors = ['#d8db20', '#197229', '#519ba3', '#240f77'] # Ugly Yellow, Green, Slate Blue, Dark Blue

colors = ['#1f77b4', '#ff7f0e', '#be52e5', '#2ca02c'] # Blue, Orange, Purple, Green

def make_plot(show = False):
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Plot
    x = frame_range
    y = np.array(extents)

    kernel = 5
    smooth_y = util.smooth(y, kernel)

    plot.plot(x, y, c = colors[1], linewidth = linewidth, alpha = alpha)
    plot.plot(x, smooth_y, c = colors[1], linewidth = linewidth)

    # Axes
    plot.xlim(x[0], x[-1])

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)
    plot.ylim(0, 360)

    # Annotate Axes
    plot.xlabel(r"$t \mathrm{\ (planet\ orbits)}$", fontsize = fontsize + 2)
    plot.ylabel(r"Azimuthal Extents $\mathrm{(degrees)}$", fontsize = fontsize + 2)

    #plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 1.0)) # outside of plot
    #plot.legend(loc = "upper left") # outside of plot

    # Title
    #title = r"$\mathrm{Azimuthal\ Extents}$"
    title = 
    plot.title("%s" % (title), y = 1.01, fontsize = fontsize + 3)

    # Save, Show, and Close
    current_directory = os.getcwd().split("/")[-3]
    current_beam = os.getcwd().split("/")[-1]
    if version is None:
        save_fn = "%s/extents-%s-%s.png" % (save_directory, current_directory, current_beam)
    else:
        save_fn = "%s/v%04d_extents-%s-%s.png" % (save_directory, version, current_directory, current_beam)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

    # Save Data
    ### TBD ###


##### Make Plots! #####

make_plot(show = show)

