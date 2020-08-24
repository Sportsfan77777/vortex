"""
plots ALMA azimuthal intensity profiles

Usage: python plotRealAlmaAzimuthlProfiles.py [name]
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

#from astropy.io import fits

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
    parser.add_argument('id_number', type = int,
                         help = 'id number of imaged system')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "almaAzimuthalProfiles",
                         help = 'save directory (default: almaAzimuthalProfiles)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", nargs = '+', type = float, default = None,
                         help = 'max_y for each frame, or same for all (default: None)')
    parser.add_argument('--profiles', dest = "num_profiles", type = int, default = 5,
                         help = 'number of profiles (default: 5)')
    parser.add_argument('-s', dest = "num_scale_heights", type = float, default = 16.0,
                         help = 'number of scale heights (default: 8.0)')

    parser.add_argument('-n', dest = "normalize", action = 'store_true', default = False,
                         help = 'normalize by max (default: normalize)')

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

# File Number
id_number = args.id_number

# Get Data
default_intensity = pickle.load(open(glob.glob("fits%03d*.p" % id_number)[0], "rb"))
deprojected_intensity = pickle.load(open(glob.glob("deprojected_fits%03d*.p" % id_number)[0]))
header = pickle.load(open(glob.glob("params_fits%03d*.p" % id_number)[0]))

rad = header["rad"]
theta = header["theta"]

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
normalize = args.normalize

show = args.show
max_y = args.max_y
if normalize:
    max_y = 1

num_profiles = args.num_profiles
num_scale_heights = args.num_scale_heights

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

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

ns = [num_scale_heights / 2, num_scale_heights / 4, num_scale_heights / 4, num_scale_heights / 2]
labels = [r"$\mathrm{-%.01f\ h}$" % ns[0], r"$\mathrm{-%.01f\ h}$" % ns[1], r"$\mathrm{+0\ h}$", r"$\mathrm{+%0.1f\ h}$" % ns[-2], r"$\mathrm{+%0.1f\ h}$" % ns[-1]]

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    px_scale = header['cdelt2'] * 3600
    num_x = header['naxis1']; num_y = header['naxis2']
    width = num_x * px_scale; height = num_y * px_scale

    xs = np.linspace(-width / 2, width / 2, num_x)
    ys = np.linspace(-height / 2, height / 2, num_x)

    ### Data ###
    intensity_cart = np.copy(deprojected_intensity)
    rs, thetas, rs_grid, thetas_grid, intensity_polar = sq.cartesian_to_polar(intensity_cart, xs, ys)
    if normalize:
        intensity_polar /= np.max(intensity_polar)

    intensity_polar = np.roll(intensity_polar, len(thetas) / 2, axis = 0) # Note: Must transpose after!!!
    azimuthal_radii, azimuthal_profiles = az.get_profiles(intensity_polar.T, header, args, shift = None, start = 0.2, end = 1)

    ### Plot ###
    # Profiles
    x = theta * (180.0 / np.pi) - 180.0
    for i, (radius, azimuthal_profile) in enumerate(zip(azimuthal_radii, azimuthal_profiles)):
        plot.plot(x, azimuthal_profile, linewidth = linewidth, c = colors[i], alpha = alpha, label = labels[i])

    # Axes
    max_x = 180
    plot.xlim(-max_x, max_x)
    angles = np.linspace(-max_x, max_x, 7)
    plot.xticks(angles)

    if max_y is None:
        plot.ylim(0, plot.ylim()[-1]) # No Input
    else:
        plot.ylim(0, max_y) # Input

    # Annotate Axes
    plot.xlabel(r"$\phi - \phi_\mathrm{0}$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)
    plot.ylabel(r"$\mathrm{Intensity}$", fontsize = fontsize)

    # Legend
    plot.legend(loc = "upper right", bbox_to_anchor = (1.34, 0.94)) # outside of plot

    # Extra Annotation
    rc_line = r"$r_\mathrm{c} = %.02f$" % azimuthal_radii[(num_profiles - 1) / 2]
    plot.text(-170, 0.90 * plot.ylim()[-1], rc_line, fontsize = fontsize, horizontalalignment = 'left')
 
    center_x = 1.38 * plot.xlim()[-1]
    top_y = plot.ylim()[-1]

    line = "Radii"
    plot.text(center_x, 0.95 * top_y, line, fontsize = fontsize, horizontalalignment = 'center')
    plot.text(center_x, 0.95 * top_y, line, fontsize = fontsize, horizontalalignment = 'center')

    # Title
    #title = "\n" + r"$t$ $=$ $%.1f$   " % (orbit) + "[$m_p(t)$ $=$ $%.2f$ $M_J$]" % (current_mass)
    #plot.title("%s" % (title), fontsize = fontsize + 1)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/almaAzimuthalProfiles_%s.png" % (save_directory, header['savename'])
    else:
        save_fn = "%s/v%04d_almaAzimuthalProfiles_%s.png" % (save_directory, version, header['savename'])
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)

