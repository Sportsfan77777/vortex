"""
plot 3-D density slices of 3-D simulations

Options:

Usage:
python plotDensityFluctuationSlices.py --frame 500 -z 0 (or -t, -r)
python plotDensityFluctuationSlices.py --frame 500 -z 0 (or -t, -r) --a 1.0 -b 2.0 # center on vortex
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

save_directory = "gasDensityFluctuationSlices"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
pluto_par = pickle.load(open(param_fn, "rb"))

num_rad = int((pluto_par["X1-grid"])[2])
num_z = int((pluto_par["X2-grid"])[2])
num_theta = int((pluto_par["X3-grid"])[2])

if ((pluto_par["X1-grid"])[3] == "l+"):
    # log grid
    rad = np.linspace(np.log10(float((pluto_par["X1-grid"])[1])), np.log10(float((pluto_par["X1-grid"])[4])), num_rad)
else:
    # uniform grid
    rad = np.linspace(float((pluto_par["X1-grid"])[1]), float((pluto_par["X1-grid"])[4]), num_rad)
zs = np.linspace(float((pluto_par["X2-grid"])[1]), float((pluto_par["X2-grid"])[4]), num_z)
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float((pluto_par["P_Sigma0"])[0])
scale_height = float((pluto_par["P_AspectRatio"])[0])

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
    normalization = surface_density_zero / (np.sqrt(2.0 * np.pi) * scale_height)

    density = (fromfile("rho.%04d.dbl" % frame).reshape(num_theta, num_z, num_rad))
    normalized_density = density / normalization

    initial_density = (fromfile("rho.0000.dbl").reshape(num_theta, num_z, num_rad))
    normalized_density_initial = initial_density / normalization

    ### Plot ###
    if o.r_slice is not None:
        # Axes
        xs = theta
        ys = zs
        # Slice
        slice_choice = "rad"; this_slice = o.r_slice
        this_slice_i = np.searchsorted(rad, this_slice)
        density_slice = normalized_density[:, :, this_slice_i].T
        initial_slice = normalized_density_initial[:, :, this_slice_i].T
    elif o.t_slice is not None:
        # Axes
        xs = rad
        ys = zs
        # Slice
        slice_choice = "phi"; this_slice = o.t_slice
        this_slice_i = np.searchsorted(theta, this_slice)
        density_slice = normalized_density[this_slice_i, :, :]
        initial_slice = normalized_density_initial[this_slice_i, :, :]
    elif o.z_slice is not None:
        # Axes
        xs = rad
        ys = theta
        # Slice
        slice_choice = "theta"; this_slice = o.z_slice
        this_slice_i = np.searchsorted(zs, this_slice)
        density_slice = normalized_density[:, this_slice_i, :]
        initial_slice = normalized_density_initial[:, this_slice_i, :]

    noise = 10**(-10)
    density_fluctutation_slice = np.log10(abs(density_slice - initial_slice) / initial_slice + noise) # Fluctuation

    result = ax.pcolormesh(xs, ys, density_fluctutation_slice, cmap = cmap) # log10 of |fluctuation|

    # Diagnostic
    print np.max(density_fluctutation_slice)

    # Colorbar
    clim = [o.clim_a, o.clim_b]

    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    # Limits
    if o.r_slice is not None:
        plot.xlim(o.t_in, o.t_out)
        plot.ylim(o.z_in, o.z_out)
    elif o.t_slice is not None:
        plot.xlim(o.r_in, o.r_out)
        plot.ylim(o.z_in, o.z_out)
    elif o.z_slice is not None:
        plot.xlim(o.r_in, o.r_out)
        plot.ylim(o.t_in, o.t_out)

    # Annotate
    r_label = "Radius"; phi_label = "Phi"; z_label = "Theta"

    if o.r_slice is not None:
        xlabel = phi_label; ylabel = z_label; suffix = "tz"
    elif o.t_slice is not None:
        xlabel = r_label; ylabel = z_label; suffix = "rz"
    elif o.z_slice is not None:
        xlabel = r_label; ylabel = phi_label; suffix = "rt"

    #this_title = readTitle()
    plot.xlabel(xlabel, fontsize = fontsize)
    plot.ylabel(ylabel, fontsize = fontsize)
    plot.title("Gas Density Fluctuation Slice (%s = %s) at Orbit %d" % (slice_choice, this_slice, orbit), fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("%s/%s_densityFluctuationSlice_%s%04d.png" % (save_directory, o.prefix, suffix, frame), bbox_inches = 'tight', dpi = my_dpi)
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
                    dest="s_start", type = "int", default = 0,
                    help="start of sample (range call)")
  parser.add_option("-w", 
                    dest="s_end", type = "int", default = max_frame,
                    help="end of sample (range call)")
  parser.add_option("-x", 
                    dest="s_rate", type = "int", default = 20,
                    help="rate of sample (range call)")

  # Which direction to skip? (Store the slice)
  parser.add_option("-r",
                    dest="r_slice", type = "float", default = None,
                    help="if missing radial direction, store r slice")
  parser.add_option("-t",
                    dest="t_slice", type = "float", default = None,
                    help="if missing azimuthal direction, store t slice")
  parser.add_option("-z",
                    dest="z_slice", type = "float", default = None,
                    help="if missing theta direction (out of the plane!), store z slice")

  # Ranges in plot? (defaults are domain ranges)
  parser.add_option("-a", 
                    dest="r_in", type = "float", default = rad[0],
                    help="start of r range")
  parser.add_option("-b", 
                    dest="r_out", type = "float", default = rad[-1],
                    help="end of r range")
  parser.add_option("-c", 
                    dest="t_in", type = "float", default = theta[0],
                    help="start of phi range")
  parser.add_option("-d", 
                    dest="t_out", type = "float", default = theta[-1],
                    help="end of phi range")
  parser.add_option("-e", 
                    dest="z_in", type = "float", default = zs[0],
                    help="start of theta range (out of the plane!)")
  parser.add_option("-f", 
                    dest="z_out", type = "int", default = zs[-1],
                    help="end of theta range (out of the plane!)")

  # Savename
  parser.add_option("--name", 
                    dest="prefix", default = "default",
                    help="savename to identify input parameters for plot")

  # Colorbar Limits
  parser.add_option("-p", 
                    dest="clim_a", type = "float", default = -2.0,
                    help="clim min")
  parser.add_option("-q", 
                    dest="clim_b", type = "float", default = 0.0,
                    help="clim max")

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



