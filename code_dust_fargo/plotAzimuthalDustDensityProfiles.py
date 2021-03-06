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
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalDensity",
                         help = 'save directory (default: azimuthalDensity)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'radial range in plot (default: None)')
    parser.add_argument('--profiles', dest = "num_profiles", type = int, default = 5,
                         help = 'number of profiles (default: 5)')
    parser.add_argument('-s', dest = "num_scale_heights", type = float, default = 1.0,
                         help = 'number of scale heights (default: 1.0)')

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
max_y = args.max_y

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

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

colors = ["#ff7f0e", "#9467bd", "#1f77b4", "#2ca02c", "#d62728"]
dashes = [[3, 3], [42, 4], [10000, 1], [42, 4], [3, 3]]

def make_plot(frame, shift, azimuthal_radii, azimuthal_profiles, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Plot ###
    x = theta * (180.0 / np.pi) - 180.0
    for i, (radius, azimuthal_profile) in enumerate(zip(azimuthal_radii, azimuthal_profiles)):
        plot.plot(x, azimuthal_profile, linewidth = linewidth, c = colors[i], dashes = dashes[i], alpha = alpha, label = "%.3f" % radius)

    # Mark Planet
    if shift is None:
        planet_loc = theta[0]
    else:
        if shift < -len(theta):
            shift += len(theta)
        planet_loc = theta[shift] * (180.0 / np.pi) - 180.0
    plot.scatter(planet_loc, 0, c = "k", s = 150, marker = "D") # planet

    # Axes
    plot.xlim(-180, 180)

    angles = np.linspace(-180, 180, 7)
    plot.xticks(angles)

    if max_y is not None:
        plot.ylim(0, max_y)
    else:
        plot.ylim(0, plot.ylim()[-1])

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    plot.xlabel(r"$\phi - \phi_\mathrm{center}$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)
    #plot.ylabel("Azimuthal Dust Density", fontsize = fontsize)
    plot.ylabel(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{dust}$", fontsize = fontsize)

    size_label = util.get_size_label(size)
    stokes_number = util.get_stokes_number(size)

    title = r"%s$\mathrm{-size}$ $\mathrm{(St}_\mathrm{0}$ $=$ $%.03f \mathrm{)}$" % (size_label, stokes_number)
    #title = r"(t = %.1f orbits)" % (orbit)
    plot.title("%s" % (title), y = 1.01, fontsize = fontsize + 1)

    # Annotate Parameters
    line_x = -165.0; line_y = 0.91; linebreak = 0.07
    center_x = 1.30 * plot.xlim()[-1]; top_y = plot.ylim()[-1]

    #left1 = r"$M_p = %d$ $M_{Jup}$" % (planet_mass)
    #left2 = r"$\nu = 10^{%d}$" % (np.log10(viscosity))
    #plot.text(line_x, line_y * plot.ylim()[-1], left1, horizontalalignment = 'left', fontsize = fontsize)
    #plot.text(line_x, (line_y - linebreak) * plot.ylim()[-1], left2, horizontalalignment = 'left', fontsize = fontsize)

    #right1 = r"$T_{growth} = %d$ $\rm{orbits}$" % (taper)
    #right2 = r"$s$ = %s" % (size_label)
    #plot.text(360 - line_x, line_y * plot.ylim()[-1], right1, horizontalalignment = 'right', fontsize = fontsize)
    #plot.text(360 - line_x, (line_y - linebreak) * plot.ylim()[-1], right2, horizontalalignment = 'right', fontsize = fontsize)

    left1 = r"$t = %.1f$" % (orbit)
    plot.text(line_x, line_y * plot.ylim()[-1], left1, horizontalalignment = 'left', fontsize = fontsize)


    plot.text(center_x, 0.95 * top_y, "Radii", fontsize = fontsize - 1, horizontalalignment = 'center')

    plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 0.94)) # outside of plot

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/azimuthalDensityProfiles_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_azimuthalDensityProfiles_%04d.png" % (save_directory, version, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

###############################################################################

def full_procedure(frame, show = False):
    """ Every Step """
    # Read Data
    density = util.read_dust_data(frame, fargo_par)

    # Choose shift option
    if center:
        """
        # Choose source
        if taper < 10.1:
            src_density = np.copy(density)
            this_threshold = threshold
        elif taper > 999.9:
            # hum-size or larger (for now, T = 1000)
            src_density = np.copy(density)
            this_threshold = threshold
        else:
            # smaller than hmm-size (right now micron only)
            hmm_directory = "/rsgrps/kkratterstudents/mhammer/fargo_tests/fargoDUST/one_jupiter_old/half-mm-size/no_diffusion/taper%d" % taper
            src_density = util.read_dust_data(frame, fargo_par, directory = hmm_directory)

            ##### Set threshold here!!!! ####
            this_threshold = util.get_threshold(0.01)
        """

        # Center vortex
        if fargo_par["MassTaper"] < 10.1:
            shift_c = az.get_azimuthal_peak(density, fargo_par)
        else:
            shift_c = az.get_azimuthal_center(density, fargo_par, threshold = threshold)
    else:
        shift_c = None

    # Get and plot profiles
    azimuthal_radii, azimuthal_profiles = az.get_profiles(density, fargo_par, args, shift = shift_c)
    make_plot(frame, shift_c, azimuthal_radii, azimuthal_profiles, show = show)

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

