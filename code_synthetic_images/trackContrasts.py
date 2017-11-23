"""
tracks vortex centers
(intended use: For T = 1000, compare each center to 'hcm-size center')
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams

import util
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])


size_names = ["cm", "hcm", "mm", "hmm", "hum", "um"]
sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot azimuthal density profiles."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = 3,
                         help = 'select range(start, end, rate)')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = ".",
                         help = 'save directory (default: current directory)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'radial range in plot (default: None)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
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
fn = "id%04d_par.p" % args.id_number
fargo_par = pickle.load(open(fn, "rb"))

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

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

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth

alpha = args.alpha
dpi = args.dpi

###############################################################################

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

fargo_par["new_num_rad"] = new_num_rad; fargo_par["new_num_theta"] = new_num_theta
fargo_par["new_r_min"] = new_r_min; fargo_par["new_r_max"] = new_r_max

###############################################################################

def get_frame_index(frame):
    """ return array index corresponding to frame """
    return np.searchsorted(frame_range, frame)

def get_contrast(intensity):
    """ return contrast """
    contrast, _, _ = az.get_contrast(intensity, fargo_par) # returns contrast, max, opposite
    
    return contrast

def save_contrasts():
    """ save array of contrasts """
    contrast_array = np.array(contrasts)

    save_fn = "%s/id%04d_contrasts_lambda%04d_beam%03d.p" % (save_directory, id_number, wavelength, beam_size)
    pickle.dump(contrast_array, open(save_fn, "wb"))

###############################################################################

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Line Plot ###

    x = frame_range
    y = np.array(contrasts)
    plot.plot(x, y, linewidth = linewidth)

    # Annotate Axes
    plot.xlabel("Time (planet orbits)", fontsize = fontsize)
    plot.ylabel("Contrasts", fontsize = fontsize)
    #plot.title("")

    # Axes
    plot.xlim(frame_range[0], frame_range[-1])
    if max_y is not None:
        plot.ylim(0, max_y)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/id%04d_contrasts_lambda%04d_beam%03d.png" % (save_directory, id_number, wavelength, beam_size)
    else:
        save_fn = "%s/v%04d_id%04d_contrasts_lambda%04d_beam%03d.png" % (save_directory, version, id_number, wavelength, beam_size)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


###############################################################################

def start():
    """ Initialize Arrays """
    global contrasts
    contrasts = mp_array('d', len(frame_range))

def full_procedure(frame):
    """ Every Step """

    intensity = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number)
    contrast_i = get_contrast(intensity)

    frame_i = get_frame_index(frame)
    contrasts[frame_i] = contrast_i

##### Gather Data! #####

# Iterate through frames

start()

if num_cores > 1:
    p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
    p.map(full_procedure, frame_range)
    p.terminate()
else:
    for frame in frame_range:
        full_procedure(frame)

##### Make Plots! #####

make_plot(show = show)
