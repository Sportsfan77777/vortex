"""
plot 2-D density slices of 3-D simulations

Options:

Usage:
python plotDensitySlices.py -r (-t, -z)
python plotDensitySlices.py
"""

import sys
import os
import subprocess
import pickle
import glob
import time
from multiprocessing import Pool

from optparse import OptionParser

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
from readTitle import readTitle

save_directory = "gasDensitySlices"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
pluto_par = pickle.load(open(param_fn, "rb"))

num_rad = int((pluto_par["X1-grid"])[2])
num_phi = int((pluto_par["X2-grid"])[2])
num_theta = int((pluto_par["X3-grid"])[2])

rad = np.linspace(float((pluto_par["X1-grid"])[1]), float((pluto_par["X1-grid"])[4]), num_rad)
phi = np.linspace(float((pluto_par["X2-grid"])[1]), float((pluto_par["X2-grid"])[4]), num_phi)
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float((pluto_par["Sigma0_Param"])[0])
scale_height = float((pluto_par["AspectRatio_Param"])[0])

max_frame = 100

### Helper Functions ###

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "RdYlBu_r"
clim = [0, 2]

fontsize = 14
my_dpi = 100


def make_plot(frame, show = False):
    # Orbit Number
    #time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    #orbit = int(round(time / (2 * np.pi), 0)) * i
    orbit = frame

    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    ax = fig.add_subplot(111)

    # Data
    density = (fromfile("rho.%04d.dbl" % frame).reshape(num_theta, num_phi, num_rad))
    normalized_density = density / surface_density_zero

    ### Plot ###
    if skip_r:
        # Axes
        xs = phi
        ys = theta
        # Slice
        this_slice = np.searchsorted(rad, skip_r)
        density_slice = normalized_density[:, :, this_slice]
    elif skip_t:
        # Axes
        xs = rad
        ys = theta
        # Slice
        this_slice = np.searchsorted(phi, skip_t)
        density_slice = normalized_density[:, this_slice, :]
    elif skip_z:
        # Axes
        xs = rad
        ys = phi
        # Slice
        this_slice = np.searchsorted(theta, skip_z)
        density_slice = normalized_density[this_slice, :, :, ]

    result = ax.pcolormesh(xs, ys, density_slice, cmap = cmap)
    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    # Limits
    if skip_r:
        plot.xlim(o.t_in, o.t_out)
        plot.ylim(o.z_in, o.z_out)
    elif skip_t:
        plot.xlim(o.r_in, o.r_out)
        plot.ylim(o.z_in, o.z_out)
    elif skip_z:
        plot.xlim(o.r_in, o.r_out)
        plot.ylim(o.t_in, o.t_out)

    x = rad
    prefix = ""
    plot.xlim(rad[0], rad[-1])
    xlabel = "Radius"

    # Annotate
    rad_label = "Radius"; phi_label = r"$\phi$"; z_label = r"$\theta$"

    if skip_r:
        xlabel = phi_label; ylabel = theta_label
    elif skip_t:
        xlabel = r_label; ylabel = z_label
    elif skip_z:
        xlabel = r_label; ylabel = theta_label

    #this_title = readTitle()
    plot.xlabel(xlabel, fontsize = fontsize)
    plot.ylabel(ylabel, fontsize = fontsize)
    plot.title("Gas Density Map at Orbit %d" % (orbit), fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("%s/%sdensityMap_%04d.png" % (save_directory, prefix, i), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

### Option Parser ###

def new_option_parser():
  parser = OptionParser()
  # Frame(s)
  parser.add_option("--frame", default = None,
                    dest="frame", type = "int",
                    help="frame to plot")
  parser.add_option("--sample", action = "store_true", default = False,
                    dest="sample",
                    help="range of frames to plot (do not show) --- use vwx to specify")
  parser.add_option("-v", 
                    dest="s_start", type = "int",
                    help="start of sample (range call)")
  parser.add_option("-w", 
                    dest="s_end", type = "int", default = max_frame,
                    help="end of sample (range call)")
  parser.add_option("-x", default = 20,
                    dest="s_rate", type = "int",
                    help="rate of sample (range call)")


  # Which direction to skip? (Store the slice)
  parser.add_option("-r",
                    dest="skip_r", type = "float", default = None,
                    help="if missing radial direction, store r slice")
  parser.add_option("-t",
                    dest="skip_t", type = "float", default = None,
                    help="if missing azimuthal direction, store t slice")
  parser.add_option("-z",
                    dest="skip_z", type = "float", default = None,
                    help="if missing theta direction (out of the plane!), store z slice")

  # Ranges in plot? (defaults are domain ranges)
  parser.add_option("-a", 
                    dest="r_in", type = "float", default = rad[0],
                    help="start of r range")
  parser.add_option("-b", 
                    dest="r_out", type = "float", default = rad[-1],
                    help="end of r range")
  parser.add_option("-c", 
                    dest="t_out", type = "float", default = phi[0],
                    help="start of phi range")
  parser.add_option("-d", 
                    dest="t_out", type = "float", default = phi[-1],
                    help="end of phi range")
  parser.add_option("-e", 
                    dest="z_in", type = "float", default = theta[0],
                    help="start of theta range (out of the plane!)")
  parser.add_option("-f", 
                    dest="z_out", type = "int", default = theta[-1],
                    help="end of theta range (out of the plane!)")

  return parser

##### Plot One File #####

o, arguments = new_option_parser().parse_args()

if o.sample:
    sample = range(o.s_start, o.s_end, o.s_rate)
    for frame_i in sample:
        make_plot(frame_i)
elif o.frame is not None:
    make_plot(o.frame, show = True)
else:
    print "Must specify frame(s)!"



