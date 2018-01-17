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
    parser.add_argument('frames', type = int, nargs = 2,
                         help = 'select two frames to compare dust density (first is for T = 10, second is for T = 1000)')

    # Directory Selection
    parser.add_argument('--dir1', dest = "directory1", default = '../taper10',
                         help = 'select first directory to compare intensity (first is for T = 10, second is for T = 10) (default: ../taper10)')
    parser.add_argument('--dir2', dest = "directory2", default = '../taper1000',
                         help = 'select second directory to compare intensity (first is for T = 10, second is for T = 1000) (default: ../taper1000)')
    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalDensityEvolution",
                         help = 'save directory (default: azimuthalDensityEvolution)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", nargs = '+', type = float, default = None,
                         help = 'max_y for each frame, or same for all (default: None)')
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

    # Analytic Profile Parameters
    parser.add_argument('--diff', dest = "diffusion_factor", type = float, default = 1,
                         help = 'extra diffusion (default: 1)')
    parser.add_argument('-r', dest = "r_a", type = float, default = None,
                         help = 'aspect ratio r (default: 1.7 [10], 1.5 [1000])')
    parser.add_argument('--dr', dest = "dr_a", type = float, default = None,
                         help = 'aspect ratio dr (default: 0.30 [10], 0.25 [1000])')
    parser.add_argument('--dtheta', dest = "dtheta_a", type = float, default = None,
                         help = 'aspect ratio d/theta (default: 120 [10], 240 [1000])')

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
frame_range = args.frames

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show
max_y = args.max_y
if max_y is None:
    pass
elif len(max_y) == 1:
    max_y = [max_y, max_y, max_y, max_y]

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

### Analytic Parameters ###
diffusion_factor = args.diffusion_factor

# Vortex radius
if args.r_a is None:
    if taper_time < 10.1:
        r_a = 1.7
    else:
        r_a = 1.5
else:
    r_a = args.r_a
# Vortex radial width
if args.dr_a is None:
    if taper_time < 10.1:
        dr_a = 0.3
    else:
        dr_a = 0.25
else:
    dr_a = args.dr_a
# Vortex azimuthal width
if args.dtheta_a is None:
    if taper_time < 10.1:
        dtheta_a = 120.0
    else:
        dtheta_a = 240.0
else:
    dtheta_a = args.dtheta_a

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

labels = [r"$\mathrm{-0.50\ h}$", r"$\mathrm{-0.25\ h}$", r"$\mathrm{+0\ h}$", r"$\mathrm{+0.25\ h}$", r"$\mathrm{+0.50\ h}$"]

def add_to_plot(frame, fig, ax, size_name, num_frames, frame_i):
    # Convert size to number
    size_name = "cm"
    size = util.get_size(size_name)

    ### Data ###
    density1 = util.read_data(frame, 'dust', fargo_par, directory = directory = "taper10/synthetic/lambda%04d/beam%03d" % (args.wavelength, args.beam_size)) / surface_density_zero
    density2 = util.read_data(frame, 'dust', fargo_par, directory = directory = "taper1000/synthetic/lambda%04d/beam%03d" % (args.wavelength, args.beam_size)) / surface_density_zero

    # Choose shift option
    if center:
        shift1 = az.get_azimuthal_peak(density1, fargo_par)

        threshold = util.get_threshold(size)
        shift2 = az.get_azimuthal_center(density2, fargo_par, threshold = threshold)
    else:
        shift1 = None; shift2 = None

    azimuthal_radii1, azimuthal_profiles1 = az.get_profiles(density1, fargo_par, args, shift = shift1)
    azimuthal_radii1, azimuthal_profiles2 = az.get_profiles(density2, fargo_par, args, shift = shift2)

    ### Plot ###
    # Profiles
    x = theta * (180.0 / np.pi) - 180.0
    for i, (radius, azimuthal_profile) in enumerate(zip(azimuthal_radii1, azimuthal_profiles1)):
        plot.plot(x, azimuthal_profile, linewidth = linewidth, c = colors[i], alpha = alpha, label = labels[i])
    for i, (radius, azimuthal_profile) in enumerate(zip(azimuthal_radii2, azimuthal_profiles2)):
        plot.plot(x, azimuthal_profile, linewidth = linewidth, c = colors[i], alpha = alpha)

    # Mark Planet
    if shift1 is None:
        planet_loc = theta[0]
    else:
        if shift < -len(theta):
            shift += len(theta)
        planet_loc = theta[shift1] * (180.0 / np.pi) - 180.0

    if shift2 is None:
        planet_loc = theta[0]
    else:
        if shift2 < -len(theta):
            shift2 += len(theta)
        planet_loc = theta[shift2] * (180.0 / np.pi) - 180.0
    plot.scatter(planet_loc1, 0, c = "r", s = 150, marker = "D", zorder = 100) # planet
    plot.scatter(planet_loc2, 0, c = "b", s = 150, marker = "D", zorder = 100) # planet

    # Axes
    max_x = 180
    plot.xlim(-max_x, max_x)
    angles = np.linspace(-max_x, max_x, 7)
    plot.xticks(angles)

    if max_y is None:
        plot.ylim(0, plot.ylim()[-1]) # No Input
    else:
        plot.ylim(0, max_y[frame_i - 1]) # Input

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    current_mass = util.get_current_mass(orbit, taper_time, planet_mass = planet_mass)

    plot.xlabel(r"$\phi - \phi_\mathrm{center}$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)

    if frame_i == 1:
        plot.ylabel(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{dust}$", fontsize = fontsize)

    # Legend
    if frame_i == 2:
        plot.legend(loc = "upper right", bbox_to_anchor = (1.34, 0.94)) # outside of plot

    # Extra Annotation
    rc_line = r"$r_\mathrm{c} = %.02f$" % azimuthal_radii[(num_profiles - 1) / 2]
    plot.text(-170, 0.90 * plot.ylim()[-1], rc_line, fontsize = fontsize, horizontalalignment = 'left')

    if frame_i == 2:    
        center_x = 1.32 * plot.xlim()[-1]
        top_y = plot.ylim()[-1]

        line1 = "Radii"
        line2 = "Analytic"
        plot.text(center_x, 0.95 * top_y, line1, fontsize = fontsize, horizontalalignment = 'center')
        plot.text(center_x, 0.25 * top_y, line2, fontsize = fontsize, horizontalalignment = 'center')

        half_width = 0.12 * plot.xlim()[-1]
        analytic_legend_y0 = 0.18 * top_y

        analytic_legend_x = [1.005 * center_x - half_width, 1.005 * center_x + half_width]
        analytic_legend_y = [analytic_legend_y0, analytic_legend_y0]
        plot.plot(analytic_legend_x, analytic_legend_y, linewidth = linewidth, c = 'k', linestyle = "--", clip_on = False)

    # Title
    title = "\n" + r"$t$ $=$ $%.1f$   " % (orbit) + "[$m_p(t)$ $=$ $%.2f$ $M_J$]" % (current_mass)
    plot.title("%s" % (title), fontsize = fontsize + 1)
    
def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (12, 5), dpi = dpi)
    gs = gridspec.GridSpec(1, 2)

    frame_str = ""
    for i, frame_i in enumerate(frame_range):
        ax = fig.add_subplot(gs[i])
        ax = add_to_plot(frame_i, fig, ax, "cm", len(frame_range), i + 1)
        ax = add_to_plot(frame_i, fig, ax, "mm", len(frame_range), i + 1)
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    #### Finish Plot ####

    size = 1.0 # cm-size only (in the future, make this a parameter?)
    size_label = util.get_size_label(size)
    stokes_number = util.get_stokes_number(size)

    title = r"$\mathrm{1\ cm-size}$ $\mathrm{(St}_\mathrm{0}$ $=$ $%.03f \mathrm{)}$" % (stokes_number)
    fig.suptitle(title, y = 0.97, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)

    # Save and Close
    plot.tight_layout()

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/azimuthalDensityEvolution_%s.png" % (save_directory, frame_str)
    else:
        save_fn = "%s/v%04d_azimuthalDensityEvolution_%s.png" % (save_directory, version, frame_str)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)
