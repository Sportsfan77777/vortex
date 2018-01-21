"""
plots ALMA fits files

Usage: python plotRealAlmaData.py [name]
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

from astropy.io import fits

import math
import numpy as np
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline as spline

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pylab import rcParams, fromfile

import util
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])


def new_argument_parser(description = "Plot real ALMA images."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('name',
                         help = 'name of imaged system')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "almaImages",
                         help = 'save directory (default: almaImages)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')


    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: viridis)')
    parser.add_argument('--cmax', dest = "cmax", type = float, default = None,
                         help = 'maximum density in colorbar (default: 10 for hcm+, 2.5 otherwise)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser


###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

name = 'SY_Cha'
filename = 'J10563044_centered.fits' # replace with args.name eventually

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show
version = args.version

# Plot Parameters (constant)
cmap = args.cmap
cmax = args.cmax
fontsize = args.fontsize
dpi = args.dpi

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    #alma_image = pyfits.getdata(filename, 0, header = True)
    fits_file = fits.open(filename)[0]
    intensity = fits_file.data[0, 0, :, :]

    #pixscales_alma = alma_image[1]['CDELT2'] * 3600

    ### Plot ###
    x = np.linspace(-10, 10, 1024)
    y = np.copy(x)
    result = ax.pcolormesh(x, y, intensity, cmap = cmap)
    
    #result.set_clim(clim[0], clim[1])

    # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
    if colorbar:
        # Only for last frame
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "8%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)
        if frame_i == 2:
            cbar.set_label("Normalized Intensity", fontsize = fontsize, rotation = 270, labelpad = 25)

    # Axes
    box = 2
    plot.xlim(-box, box)
    plot.ylim(-box, box)
    plot.axes().set_aspect('equal')

    # Annotate Axes
    plot.xlabel(r"$\mathrm{Relative\ R.A.\ [arcsec]}$", fontsize = fontsize)
    plot.ylabel(r"$\mathrm{Relative\ Dec.\ [arcsec]}$", fontsize = fontsize)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/almaImage_%s.png" % (save_directory, name)
    else:
        save_fn = "%s/v%04d_almaImage_%s.png" % (save_directory, version, name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)
