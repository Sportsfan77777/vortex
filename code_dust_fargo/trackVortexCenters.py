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

### Get Input Parameters ###

# Frames
start = args.frames[0]; end = args.frames[1]; rate = args.frames[2]
frame_range = range(start, end + 1, rate)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show
max_y = args.max_y

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

###############################################################################

### Helper Methods ###
def get_fargo_par(directory):
    """ return fargo parameter dictionary """
    fargo_par = util.get_pickled_parameters(directory = directory)

    ### Create new parameters ###
    rad = np.linspace(fargo_par["Rmin"], fargo_par["Rmax"], fargo_par["Nrad"])
    theta = np.linspace(0, 2 * np.pi, fargo_par["Nsec"])

    ### Add new parameters to dictionary ###
    fargo_par["rad"] = rad
    fargo_par["theta"] = theta

    return fargo_par

def get_frame_index(frame):
    """ return array index corresponding to frame """
    return np.searchsorted(frame_range, frame)

def get_center(density, fargo_par, size):
    """ return theta for vortex center (in degrees) """
    ######## Get Parameters #########
    theta = fargo_par["theta"]

    ########### Method ############## 
    threshold = util.get_threshold(size)
    shift_c = az.get_azimuthal_center(density, fargo_par, threshold = threshold)
    theta_c = theta[shift_c] * (180.0 / np.pi)
    
    return theta_c

def angle_difference(angles1, angles2):
    """ return difference between two arrays of angles """
    angle_diff = angles1 - angles2

    # Correct for angle pairs on opposite sides of zero
    angle_diff[angle_diff > 180] -= 360
    angle_diff[angle_diff < 180] += 360

    return angle_diff

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
    reference = np.array(center_arrays["hcm"])

    for i, size_name in enumerate(size_names):
        x = frame_range
        y = angle_difference(np.array(center_arrays[size_name]), reference)
        plot.plot(x, y, c = colors[i], linewidth = linewidth, label = util.get_size_label(sizes[i]))

    # Annotate Axes
    plot.xlabel("Time (planet orbits)", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)
    plot.title("Vortex Centers")

    plot.legend(loc = "upper right")

    # Axes
    plot.xlim(frame_range[0], frame_range[-1])
    if max_y is not None:
        plot.ylim(-max_y, max_y)

    # Save, Show, and Close
    save_fn = "vortexCenters.png"
    plot.savefig(save_fn, bbox_inches = 'tight')

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


###############################################################################

def start():
    """ Initialize Arrays """
    global center_arrays
    center_arrays = {}

    for i, size_name in enumerate(size_names):
        center_arrays[size_name] = mp_array('d', len(frame_range))

def full_procedure(frame):
    """ Every Step """

    for i, size_name in enumerate(size_names):
        directory = "../%s-size/" % size_name
        fargo_par_i = get_fargo_par(directory)

        util.get_pickled_parameters(directory = directory)
        density = util.read_data(frame, 'dust', fargo_par_i, directory = directory)

        frame_i = get_frame_index(frame)
        (center_arrays[size_name])[frame_i] = get_center(density, fargo_par_i, sizes[i])

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

