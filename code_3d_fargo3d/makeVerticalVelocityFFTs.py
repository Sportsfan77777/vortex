"""
frequency vs. radius plot showing modes of the VSI from the vertical velocity
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

from scipy import signal
from scipy.interpolate import interp1d as interpolate
from scipy.ndimage import map_coordinates

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
import square as sq

from advanced import Parameters
#from reader import Fields

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('--range', dest = "frames", type = int, nargs = '+', default = None,
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('--rate', dest = "rate", type = int, default = 1,
                         help = 'frame rate (default: 1)')
    parser.add_argument('--cadence', dest = "cadence", type = int, default = 1,
                         help = 'frame rate (default: 1)')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--load', dest = "load_directory", default = "midplaneVerticalVelocityMapsOverTime",
                         help = 'load directory (default: midplaneVerticalVelocityMapsOverTime)')
    parser.add_argument('--dir', dest = "save_directory", default = "verticalVelocityFFTs",
                         help = 'save directory (default: verticalVelocityFFTs)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "hot",
                         help = 'color map (default: hot)')
    parser.add_argument('--cmax', dest = "cmax", type = float, default = 0.1,
                         help = 'min and max values in colorbar (default: 0.1)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 17,
                         help = 'fontsize of plot annotations (default: 17)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 15,
                         help = 'fontsize of plot annotations (default: 15)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
p = Parameters()

num_rad = p.ny; num_theta = p.nx; num_z = p.nz
r_min = p.ymin; r_max = p.ymax
z_min = p.zmin; z_max = p.zmax

surface_density_zero = p.sigma0

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

scale_height = p.aspectratio
viscosity = 1e-7 #p.nu

#######################################################

### Get Input Parameters ###

# Frames
if args.frames is not None:
    frame_range = util.get_frame_range(args.frames)

rate = args.rate
cadence = args.cadence

# Number of Cores 
num_cores = args.num_cores

# Files
load_directory = args.load_directory
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

show = args.show
version = args.version

cmap = args.cmap
if args.cmax is None:
    clim = [0, 0.01]
else:
    clim = [0, args.cmax]

fontsize = args.fontsize
dpi = args.dpi

rc['xtick.labelsize'] = args.labelsize
rc['ytick.labelsize'] = args.labelsize

### Helper Functions ###

def process_data(data, frames, radii):
    # Set up output
    half = len(frames) / 2 + 1
    fft_data = np.zeros((half, len(radii)))

    for i, radius in enumerate(radii):
        # One radius at a time
        sliver = data[:, i]

        # FFT
        fft_sliver = np.fft.fft(sliver, axis = -1)
        fft_data[:, i] = (np.abs(fft_sliver))[:half]
    freq = np.fft.fftfreq(sliver.shape[-1])[:half]

    return fft_data, freq

# Load Data
#if args.name is None:
directory_name = os.getcwd().split("/")[-1]
name = directory_name

data = np.array(pickle.load(open("%s/%s_verticalVelocityMap-data.p" % (load_directory, name), "rb")))[::rate, :]
frames = np.array(pickle.load(open("%s/%s_verticalVelocityMap-frames.p" % (load_directory, name), "rb")))[::rate]
radii = np.array(pickle.load(open("%s/%s_verticalVelocityMap-radii.p" % (load_directory, name), "rb")))

# Process!
fft_data, freq = process_data(data, frames, radii)
freq *= cadence

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (10, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Plot ###
    x = radii
    y = freq
    result = plot.pcolormesh(x, y, np.abs(fft_data), cmap = cmap)
    result.set_clim(clim)
    plot.colorbar()

    # Reference Line
    xref = radii
    yref1 = np.power(radii, -1.5) * (np.sqrt(5) - 2.0)
    yref2 = np.power(radii, -1.5)
    ref_color = "deepskyblue"
    plot.plot(xref, yref1, c = ref_color, linewidth = 2, linestyle = "--")
    plot.plot(xref, yref2, c = ref_color, linewidth = 2)

    # Axes
    plot.xlim(radii[0], radii[-1])
    plot.ylim(min(freq[freq > 0]), max(freq))
    plot.xscale("log")
    plot.yscale("log")

    xticks = np.concatenate((np.logspace(np.log10(0.5), np.log10(1.0), 3), np.logspace(np.log10(1.0), np.log10(2.5), 5)[1:]))
    plot.minorticks_off() # Fixes known bug where xticks aren't removed if scale is log
    plot.xticks(xticks, ['%.2f' % xtick for xtick in xticks])

    yticks = np.logspace(np.log10(1.0 * cadence / len(frames)), np.log10(0.5 * cadence), 10)
    plot.yticks(yticks, ['%.3f' % ytick for ytick in yticks])

    # Annotate Axes
    radius_unit = r'$r_\mathrm{p}$'
    plot.xlabel("Radius [%s]" % radius_unit, fontsize = fontsize)
    plot.ylabel(r"$\omega / \Omega_\mathrm{0}$", fontsize = fontsize)
    
    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1]

    if version is None:
        save_fn = "%s/%s_FFTverticalVelocityMap_%04d_%04d_%04d.png" % (save_directory, directory_name, frames[0], frames[-1], args.rate)
    else:
        save_fn = "%s/v%04d_%s_FFTverticalVelocityMap_%04d_%04d_%04d.png" % (save_directory, version, directory_name, frames[0], frames[-1], args.rate)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi, pad_inches = 0.2)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)
