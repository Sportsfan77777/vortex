"""
plot 2-D density maps

python plotDensityMaps.py
python plotDensityMaps.py frame_number
python plotDensityMaps.py -1 <<<===== Plots a sample
python plotDensityMaps.py -m
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
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import gridspec

from pylab import rcParams
from pylab import fromfile

import util
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Size Selection
    parser.add_argument('--sizes', dest = "sizes", nargs = 4, default = ["um", "cm", "hcm", "mm"],
                         help = 'select 4 sizes (default: [um, cm, hcm, mm])')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalDensityTwoByTwo",
                         help = 'save directory (default: azimuthalDensityTwoByTwo)')

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
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 18,
                         help = 'fontsize of plot annotations (default: 18)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 15,
                         help = 'labelsize of plot annotations (default: 15)')
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

### Get FARGO Parameters ###
fargo_par = util.get_pickled_parameters(directory = "../cm-size")

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Sizes
sizes = args.sizes

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show
max_y = args.max_y

num_profiles = args.num_profiles
num_scale_heights = args.num_scale_heights

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version

center = args.center
threshold = args.threshold
#if threshold is None:
#    threshold = util.get_threshold(size)

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

def add_to_plot(frame, fig, ax, size_name, num_sizes, frame_i):
    # Convert size to number
    size = util.get_size(size_name)

    ### Data ###
    density = util.read_data(frame, 'dust', fargo_par, directory = "../%s-size" % size_name) / surface_density_zero

    # Choose shift option
    if center:
        # Center vortex
        if fargo_par["MassTaper"] < 10.1:
            shift = az.get_azimuthal_peak(density, fargo_par)
        else:
            threshold = util.get_threshold(size)
            shift = az.get_azimuthal_center(density, fargo_par, threshold = threshold)
    else:
        shift = None
    azimuthal_radii, azimuthal_profiles = az.get_profiles(density, fargo_par, args, shift = shift)

    ### Plot ###
    # Profiles
    x = theta * (180.0 / np.pi) - 180.0
    for i, (radius, azimuthal_profile) in enumerate(zip(azimuthal_radii, azimuthal_profiles)):
        plot.plot(x, azimuthal_profile, linewidth = linewidth, c = colors[i], alpha = alpha, label = "%.3f" % radius)

    # Analytic
    if frame_i != 1:
        middle_i = (num_profiles - 1) / 2
        radius = azimuthal_radii[middle_i] # middle
        center_density = azimuthal_profile[middle_i][(len(azimuthal_profile[middle_i]) - 1) / 2]

        aspect_ratio = (1.5 / 0.25) * (240 * np.pi / 180.0) # (r / dr) * d\theta
        S = util.get_stokes_number(size) / (viscosity / scale_height**2) # St / \alpha

        analytic = np.array([az.get_analytic_profile(angle, aspect_ratio, S, radius) for angle in x])
        analytic = analytic / np.max(analytic) * center_density # Normalize and re-scale to density at the center
        plot.plot(x, analytic, linewidth = linewidth, alpha = alpha, linestyle = "--", c = "k")

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

    if frame_i <= 2:
        # Remove unless bottom
        ax.set_xticklabels([])

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    if frame_i > 2:
        plot.xlabel(r"$\phi - \phi_\mathrm{center}$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)

    if frame_i == 1:
        plot.ylabel(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{gas}$", fontsize = fontsize)
    elif frame_i == 3:
        plot.ylabel(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{dust}$", fontsize = fontsize)

    # Legend
    if frame_i == 2:
        plot.legend(loc = "upper right", bbox_to_anchor = (1.3, 0.94)) # outside of plot

    # Extra Annotation
    if frame_i % 2 == 0:    
        center_x = 1.32 * plot.xlim()[-1]
        top_y = plot.ylim()[-1]

        if frame_i == 2:
            line = "Radii"
        elif frame_i == 4:
            line = "Analytic"
        plot.text(center_x, 0.95 * top_y, line, fontsize = fontsize, horizontalalignment = 'center')
        plot.text(center_x, 0.95 * top_y, line, fontsize = fontsize, horizontalalignment = 'center')

        if frame_i == 4:
            half_width = 0.12 * plot.xlim()[-1]
            analytic_legend_y0 = 0.85 * top_y

            analytic_legend_x = [1.005 * center_x - half_width, 1.005 * center_x + half_width]
            analytic_legend_y = [analytic_legend_y0, analytic_legend_y0]
            plot.plot(analytic_legend_x, analytic_legend_y, linewidth = linewidth, c = 'k', linestyle = "--", clip_on = False)

    # Title
    size_label = util.get_size_label(size)
    stokes_number = util.get_stokes_number(size)

    title = r"%s$\mathrm{-size}$ $\mathrm{(St}_\mathrm{0}$ $=$ $%.03f \mathrm{)}$" % (size_label, stokes_number)
    if frame_i == 1:
        title = r"$\mathrm{Gas\ Density}$"
    plot.title("%s" % (title), fontsize = fontsize + 1)

    
def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (12, 7), dpi = dpi)
    gs = gridspec.GridSpec(2, 2)

    size_str = ""
    for i, size_i in enumerate(sizes):
        ax = fig.add_subplot(gs[i])
        ax = add_to_plot(frame, fig, ax, size_i, len(frame_range), i + 1)
        size_str += "size_i_"
    size_str = size_str[:-1] # Trim last '_'

    #### Finish Plot ####

    # Title
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    current_mass = util.get_current_mass(orbit, taper_time, planet_mass = planet_mass)
    if orbit >= 1: # taper_time:
        frame_title = r"$t$ $=$ $%.1f$" % (orbit)
        fig.suptitle(frame_title, y = 0.99, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)
    else:
        frame_title = "\n" + r"$t$ $=$ $%.1f$" % (orbit) + "\n" + "[$m_p(t)$ $=$ $%.2f$ $M_J$]" % (current_mass)
        fig.suptitle(frame_title, y = 0.99, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)

    # Save and Close
    plot.tight_layout()

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/azimuthalDensityProfiles_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_azimuthalDensityProfiles_%04d.png" % (save_directory, version, frame)
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
