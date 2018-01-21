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

    return parser


filename = 'J10563044_centered.fits'

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    #alma_image = pyfits.getdata(filename, 0, header = True)
    fits_file = fits.open(filename)[0]
    intensity = fits_file.data

    #pixscales_alma = alma_image[1]['CDELT2'] * 3600

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi)
    result = ax.pcolormesh(x, y, np.transpose(intensity), cmap = cmap)

    fig.colorbar(result)
    #result.set_clim(clim[0], clim[1])

    # Axes
    #plot.xlim(x_min, x_max)
    #plot.ylim(0, 360)

    # Annotate Axes
    plot.xlabel(r"$\mathrm{Relative R.A. [arcsec]}$", fontsize = fontsize)
    plot.ylabel(r"$\mathrm{Relative Dec. [arcsec]}$", fontsize = fontsize)

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
