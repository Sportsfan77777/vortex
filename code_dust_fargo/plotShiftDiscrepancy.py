"""
plot azimuthal density profiles
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

from pylab import rcParams
from pylab import fromfile

import util
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot azimuthal density profiles."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "shiftDiscrepancies",
                         help = 'save directory (default: shiftDiscrepancies)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'radial range in plot (default: None)')
    parser.add_argument('-s', dest = "num_scale_heights", type = float, default = 4.0,
                         help = 'number of scale heights (default: 4.0)')

    parser.add_argument('--shift_off', dest = "center", action = 'store_false', default = True,
                         help = 'do not center frame on vortex peak or middle (default: shift to center)')
    parser.add_argument('-t', dest = "threshold", type = float, default = None,
                         help = 'threshold for centering vortex with its center (default: varies with size)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
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

### Get Fargo Parameters ###
fargo_par = util.get_pickled_parameters()

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]
taper = fargo_par["MassTaper"]

size = fargo_par["PSIZE"]
size_label = util.get_size_label(size)

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
max_y = args.max_y
num_scale_heights = args.num_scale_heights

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version

center = args.center
threshold = args.threshold
if threshold is None:
    threshold = util.get_threshold(size)

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

###############################################################################

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

outer_start = 1.1; outer_end = 2.5

def find_shift(density1, density2, fargo_par, start = outer_start, end = outer_end):
    """ return shift needed to shift vortex center to 180 degrees """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############

    ### Get sliver ###

    def get_sliver(density):
        ### Identify center using threshold ###
        # Search outer disk only
        outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
        outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
        density_segment = density[outer_disk_start : outer_disk_end]

        # Get peak in azimuthal profile
        avg_density = np.average(density_segment, axis = 1) # avg over theta
        segment_arg_peak = np.argmax(avg_density)
        arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
        peak_rad = rad[arg_peak]

        # Zoom in on peak
        half_width = 0.5 * num_scale_heights * scale_height
        zoom_start = np.searchsorted(rad, peak_rad - half_width)
        zoom_end = np.searchsorted(rad, peak_rad + half_width)

        density_sliver = density[zoom_start : zoom_end]
        return density_sliver

    density_sliver1 = get_sliver(density1)
    density_sliver2 = get_sliver(density2)

    ### Test shifts ###
    theta_degrees = theta * (180.0 / np.pi)

    min_shift = 0; max_shift = 359; num_shifts = 360
    possible_shifts = np.linspace(min_shift, max_shift, num_shifts)

    mass_differences = np.zeros(num_shifts)

    for i, shift_i in enumerate(possible_shifts):
        shift = np.searchsorted(theta_degrees, shift_i)
        tmp_density_sliver1 = np.roll(density_sliver1, shift)

        diff = np.abs(density_sliver2 - tmp_density_sliver1)
        mass_differences[i] = np.sum(diff)

    return mass_differences

###############################################################################

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']


def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Data ###
    min_shift = 0; max_shift = 359; num_shifts = 360
    possible_shifts = np.linspace(min_shift, max_shift, num_shifts)

    grain1 = "cm"; grain2 = "hcm"
    density1 = (fromfile("../%s-size/gasdens%d.dat" % (grain1, frame)).reshape(num_rad, num_theta))
    density2 = (fromfile("../%s-size/gasdens%d.dat" % (grain2, frame)).reshape(num_rad, num_theta))

    ### Plot ###
    x = possible_shifts
    y = mass_differences(density1, density2)
    plot.plot(x, y, linewidth = linewidth, c = "k")

    # Axes
    plot.xlim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.xticks(angles)

    if max_y is not None:
        plot.ylim(0, max_y)
    else:
        plot.ylim(0, plot.ylim()[-1])

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    plot.xlabel(r"$\phi$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)
    #plot.ylabel("Azimuthal Dust Density", fontsize = fontsize)
    plot.ylabel(r"$\Sigma_\mathrm{diff}$", fontsize = fontsize)

    #size_label = util.get_size_label(size)
    #stokes_number = util.get_stokes_number(size)

    #title = r"%s$\mathrm{-size}$ $\mathrm{(St}_\mathrm{0}$ $=$ $%.03f \mathrm{)}$" % (size_label, stokes_number)
    #title = r"(t = %.1f orbits)" % (orbit)
    #plot.title("%s" % (title), y = 1.01, fontsize = fontsize + 1)

    #plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 0.94)) # outside of plot

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/shiftDiscrepancies_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_shiftDiscrepancies_%04d.png" % (save_directory, version, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

# Iterate through frames

if len(frame_range) == 1:
    make_plot(frame_range[0], show = show)
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(make_plot, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            make_plot(frame)

