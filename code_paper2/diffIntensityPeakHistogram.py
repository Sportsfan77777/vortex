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


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # File Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('--id', dest = "id_numbers", type = int, nargs = 2, default = None,
                         help = 'two id numbers (up to 4 digits) for these two sets of input parameters (default: None)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalIntensityEvolution",
                         help = 'save directory (default: azimuthalIntensityEvolution)')

    # Print
    parser.add_argument('--print', dest = "print_data", action = 'store_true', default = False,
                         help = 'print data (default: do not print)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    
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
fn = "id%04d_par.p" % (args.id_numbers[0])
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

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

### Get Input Parameters ###

# File Selection
frame_range = util.get_frame_range(args.frames)

id1 = args.id_numbers[0]
id2 = args.id_numbers[1]

id_str = "%04d-%04d" % (id1, id2)

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Print
print_data = args.print_data

# Plot Parameters (variable)
show = args.show

version = args.version

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

###############################################################################

### Data ###

data1 = pickle.load(open("id%04d_b%02d_intensityPeaks.p" % (id1, beam_size * planet_radius), "rb"))
data2 = pickle.load(open("id%04d_b%02d_intensityPeaks.p" % (id2, beam_size * planet_radius), "rb"))

diff = data1 - data2

if print_data:
    for (frame, d1, d2, d) in zip(frame_range, data1, data2, diff):
        print "Frame %04d: %02.1f, %02.1f (Difference: %02.2f)"

###############################################################################

##### PLOTTING #####

def make_plot(show = False):
    fig = plot.figure(figsize = (7, 5), dpi = dpi)
    ax = fig.add_subplot(111)

    # Plot
    bins = np.linspace(-20, 20, 21)
    bins = np.linspace(-20, 20, 201)
    data = diff
    plot.hist(data / len(frame_range), bins = bins1, histtype = 'step')
    plot.hist(data / len(frame_range), bins = bins2, histtype = 'step', cumulative = True)

    # Save, Show, and Close
    frame_str = ""
    for i, frame_i in enumerate(args.frames):
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    if version is None:
        save_fn = "%s/diffIntensityPeakHistogram_%s.png" % (save_directory, frame_str)
    else:
        save_fn = "%s/v%04d_diffIntensityPeakHistogram_%s.png" % (save_directory, version, frame_str)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)


