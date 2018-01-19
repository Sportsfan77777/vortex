"""
plots ALMA fits files

Usage: python plotRealAlmaData.py [name]
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

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
import pyfits as fits

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




##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    density = (fromfile("gasddens%d.dat" % frame).reshape(num_rad, num_theta))
    if center:
        if taper < 10.1:
            shift_c = az.get_azimuthal_peak(density, fargo_par)
        else:
            threshold = util.get_threshold(size)
            shift_c = az.get_azimuthal_center(density, fargo_par, threshold = threshold * surface_density_zero)
        density = np.roll(density, shift_c)
    normalized_density = density / surface_density_zero

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi)
    result = ax.pcolormesh(x, y, np.transpose(normalized_density), cmap = cmap)

    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    title = readTitle()

    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)

    if title is None:
        plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    else:
        plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/dustDensityMap_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_dustDensityMap_%04d.png" % (save_directory, version, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)
