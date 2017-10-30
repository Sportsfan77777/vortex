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
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalDensity",
                         help = 'save directory (default: gasDensityMaps)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--profiles', dest = "num_profiles", type = int, default = 5,
                         help = 'number of profiles (default: 5)')
    parser.add_argument('-s', dest = "num_scale_heights", type = float, default = 0.5,
                         help = 'number of scale heights (default: 0.5)')
    
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

### Get Fargo Parameters ###
fargo_par = util.get_pickled_parameters()

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]

size = fargo_par["PSIZE"]

### Get Input Parameters ###

# Frames
if len(args.frames) == 1:
    frame_range = args.frames
elif len(args.frames) == 3:
    start = args.frames[0]; end = args.frames[1]; rate = args.frames[2]
    frame_range = range(start, end + 1, rate)
else:
    print "Error: Must supply 1 or 3 frame arguments\nWith one argument, plots single frame\nWith three arguments, plots range(start, end + 1, rate)"
    exit()

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

version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

### Helper Methods ###

def find_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start:])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_density

def find_min(averagedDensity, peak_rad):
    try:
        outer_disk_start = np.searchsorted(rad, 1.0) # look for max radial density beyond r = 1.1
        outer_disk_end = np.searchsorted(rad, peak_rad)
        min_rad_outer_index = np.argmin(averagedDensity[outer_disk_start : outer_disk_end])

        min_index = outer_disk_start + min_rad_outer_index
        min_rad = rad[min_index]
        min_density = averagedDensity[min_index]

        #print "Min", min_rad, min_density
        return min_rad, min_density
    except:
        # No Gap Yet
        return peak_rad, 0

def find_center(density, threshold_value = 0.05):
    """ return shift needed to shift vortex center to 180 degrees """
    ### Identify center using threshold ###
    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    # Get peak in azimuthal profile
    avg_density = np.average(density_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_density)
    arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
    peak_rad = rad[arg_peak]

    # Zoom in on peak --- Average over half a scale height
    half_width = 0.25 * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    density_sliver = density[zoom_start : zoom_end]
    avg_density_sliver = np.average(density_sliver, axis = 0) # avg over rad

    # Move Minimum to Zero Degrees (vortex cannot cross zero)
    arg_min = np.argmin(avg_density_sliver)
    shift_min = int(0 - arg_min)
    avg_density_sliver = np.roll(avg_density_sliver, shift_min)

    # Spot two threshold crossovers
    threshold = threshold_value * surface_density_zero

    left_edge = np.searchsorted(avg_density_sliver, threshold)
    right_edge = len(theta) - np.searchsorted(avg_density_sliver[::-1], threshold) - 1

    center = (left_edge + right_edge) / 2.0

    ### Calculate shift for true center to 180 degrees ###
    middle = np.searchsorted(theta, np.pi)
    shift_c = int(middle - (center - shift_min))

    return shift_c

### Data ###

def get_data(frame):
    """ Gather azimuthal radii and profiles """
    # Find Peak in Radial Profile (in Outer Disk)
    density = (fromfile("gasddens%d.dat" % frame).reshape(num_rad, num_theta))
    density = density / surface_density_zero

    averagedDensity = np.average(density, axis = 1)
    peak_rad, peak_density = find_peak(averagedDensity)
    min_rad, min_density = find_min(averagedDensity, peak_rad)

    # Shift to Center
    shift_c = find_center(density, threshold_value = 5)
    density = np.roll(density, shift_c)

    # Gather Azimuthal Profiles
    num_profiles = args.num_profiles
    spread = args.num_scale_heights * scale_height / 2.0

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = [density[azimuthal_index, :] for azimuthal_index in azimuthal_indices]

    return azimuthal_radii, azimuthal_profiles

###############################################################################

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

def make_plot(frame, azimuthal_radii, azimuthal_profiles, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Plot ###
    x = theta * (180.0 / np.pi)
    for i, (radius, azimuthal_profile) in enumerate(zip(azimuthal_radii, azimuthal_profiles)):
        plot.plot(x, azimuthal_profile, linewidth = linewidth, c = colors[i], alpha = alpha, label = "%.3f" % radius)

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    plot.xlabel(r"$\phi$", fontsize = fontsize + 2)
    plot.ylabel("Density", fontsize = fontsize)
    plot.title("Azimuthal Dust Density [%f cm] \n(t = %.1f)" % (size, orbit), fontsize = fontsize + 1)

    plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 1.0)) # outside of plot)

    # Axes
    plot.xlim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.xticks(angles)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/azimuthalDensityProfiles_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_azimuthalDensityProfiles_%04d.png" % (save_directory, version, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


def full_procedure(frame, show = False):
    """ Every Step """
    azimuthal_radii, azimuthal_profiles = get_data(frame)
    make_plot(frame, azimuthal_radii, azimuthal_profiles, show = show)

##### Make Plots! #####

# Iterate through frames

if len(frame_range) == 1:
    full_procedure(frame_range[0], show = show)
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(full_procedure, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            full_procedure(frame)

