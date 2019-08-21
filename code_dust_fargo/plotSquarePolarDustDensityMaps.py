"""
plot 2-D dust density maps

usage: plotDustDensityMaps.py [-h] [-c NUM_CORES] [--dir SAVE_DIRECTORY]
                              [--hide] [-v VERSION] [--range R_LIM R_LIM]
                              [--shift] [--cmap CMAP] [--cmax CMAX]
                              [--fontsize FONTSIZE] [--dpi DPI]
                              frames [frames ...]

positional arguments:
  frames                select single frame or range(start, end, rate). error
                        if nargs != 1 or 3

optional arguments:
  -h, --help            show this help message and exit
  -c NUM_CORES          number of cores (default: 1)
  --dir SAVE_DIRECTORY  save directory (default: dustDensityMaps)
  --hide                for single plot, do not display plot (default: display
                        plot)
  -v VERSION            version number (up to 4 digits) for this set of plot
                        parameters (default: None)
  --range R_LIM R_LIM   radial range in plot (default: [r_min, r_max])
  --shift               center frame on vortex peak or middle (default: do not
                        center)
  --cmap CMAP           color map (default: viridis)
  --cmax CMAX           maximum density in colorbar (default: 10 for hcm+, 2.5
                        otherwise)
  --fontsize FONTSIZE   fontsize of plot annotations (default: 16)
  --dpi DPI             dpi of plot annotations (default: 100)
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
import azimuthal as az
import square as sq
from readTitle import readTitle

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
    parser.add_argument('--dir', dest = "save_directory", default = "squarePolarDustDensityMaps",
                         help = 'save directory (default: squarePolarDustDensityMaps)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--box', dest = "box", type = float, default = 2.5,
                         help = 'width of box (in r_p) (default: 2.5)')
    parser.add_argument('--shift', dest = "center", default = "off",
                         help = 'center frame on vortex peak or middle (default: do not center)')
    parser.add_argument('--clear', dest = "clear", action = 'store_true', default = False,
                         help = 'clear inner disk or not (default: do not clear)')

    # Plot Parameters (contours)
    parser.add_argument('--contour', dest = "use_contours", action = 'store_true', default = False,
                         help = 'use contours or not (default: do not use contours)')
    parser.add_argument('--low', dest = "low_contour", type = float, default = 1.1,
                         help = 'lowest contour (default: 1.1)')
    parser.add_argument('--high', dest = "high_contour", type = float, default = 3.5,
                         help = 'highest contour (default: 3.5)')
    parser.add_argument('--num_levels', dest = "num_levels", type = int, default = None,
                         help = 'number of contours (choose this or separation) (default: None)')
    parser.add_argument('--separation', dest = "separation", type = float, default = 0.1,
                         help = 'separation between contours (choose this or num_levels) (default: 0.1)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: inferno)')
    parser.add_argument('--cmax', dest = "cmax", type = float, default = 50,
                         help = 'maximum density in colorbar (default: 50)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
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
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
try:
  viscosity = fargo_par["Viscosity"]
except:
  viscosity = fargo_par["AlphaViscosity"]

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

box = args.box
center = args.center
clear = args.clear

# Plot Parameters (contours)
use_contours = args.use_contours
low_contour = args.low_contour
high_contour = args.high_contour
num_levels = args.num_levels
if num_levels is None:
    separation = args.separation
    num_levels = int(round((high_contour - low_contour) / separation + 1, 0))

# Plot Parameters (constant)
cmap = args.cmap
clim = [0, args.cmax]

fontsize = args.fontsize
dpi = args.dpi

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Helper Functions ###

def clear_inner_disk(data):
    # get rid of inner disk (r < outer_limit)
    filtered_data = np.copy(data)

    outer_limit = np.searchsorted(rad, 1.05)
    filtered_data[:outer_limit] = 0

    return filtered_data

def shift_density(normalized_density, fargo_par, option = "away", reference_density = None, frame = None):
    """ shift density based on option """
    if reference_density is None:
       reference_density = normalized_density

    # Options
    if option == "peak":
       shift_c = az.get_azimuthal_peak(reference_density, fargo_par)
    elif option == "threshold":
       threshold = util.get_threshold(fargo_par["PSIZE"])
       shift_c = az.get_azimuthal_center(reference_density, fargo_par, threshold = threshold)
    elif option == "away":
       shift_c = az.shift_away_from_minimum(reference_density, fargo_par)
    elif option == "lookup":
       shift_c = az.get_lookup_shift(frame)
    else:
       print "Invalid centering option. Choose (cm-)peak, (cm-)threshold, (cm-)away, or lookup"

    # Shift
    shifted_density = np.roll(normalized_density, shift_c)
    return shifted_density, shift_c

###############################################################################

def generate_colors(n):
    c = ['yellow', 'b', 'firebrick', 'w', 'green']
    colors = []
    for i in range(n):
        colors.append(c[i % len(c)])
    return colors

##### PLOTTING #####

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    density = (fromfile("gasddens%d.dat" % frame).reshape(num_rad, num_theta))
    gas_density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    normalized_density = density / (surface_density_zero / 100.0)

    # Shift
    if center is "off":
       shift_c = 0
    #elif center.startswith("cm"):
    #   normalized_density, shift_c = shift_density(normalized_density, fargo_par, option = center[3:], reference_density = cm_dust_density, frame = frame)
    else:
       normalized_density, shift_c = shift_density(normalized_density, fargo_par, option = center, frame = frame)

    # Clear
    if clear:
       normalized_density = clear_inner_disk(normalized_density)

    # Convert to cartesian
    xs, ys, _, _, normalized_density_cart = sq.polar_to_cartesian(normalized_density, rad, theta)

    ### Plot ###
    result = ax.pcolormesh(xs, ys, np.transpose(normalized_density_cart), cmap = cmap)

    cbar = fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    cbar.set_label(r"Normalized Surface Density", fontsize = fontsize + 2, rotation = 270, labelpad = 20)

    if use_contours:
        levels = np.linspace(low_contour, high_contour, num_levels)
        colors = generate_colors(num_levels)
        plot.contour(x, y, np.transpose(normalized_density_cart), levels = levels, origin = 'upper', linewidths = 1, colors = colors)

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    fig.gca().add_artist(circle)

    # Label star and planet
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    planet_size = (current_mass / planet_mass)
    plot.scatter(0, 0, c = "white", s = 300, marker = "*", zorder = 100) # star
    plot.scatter(0, 1, c = "white", s = int(70 * planet_size), marker = "D") # planet

    # Add minor grid lines
    alpha = 0.2
    dashes = [1, 5]
    #plot.grid(b = True, which = 'major', color = "black", dashes = dashes, alpha = alpha)
    #plot.grid(b = True, which = 'minor', color = "black", dashes = dashes, alpha = alpha)
    #plot.minorticks_on()

    # Axes
    plot.xlim(-box, box)
    plot.ylim(-box, box)
    plot.axes().set_aspect('equal')

    # Annotate Axes
    unit = "r_\mathrm{p}"
    plot.xlabel(r"$x$ [$%s$]" % unit, fontsize = fontsize)
    plot.ylabel(r"$y$ [$%s$]" % unit, fontsize = fontsize)

    title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)
    plot.text(0.0, 3.24, title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    plot.text(-1.7, 3.24, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    plot.text(1.7, 3.24, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')
    # Save, Show, and Close
    if version is None:
        save_fn = "%s/squareDensityMap_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_squareDensityMap_%04d.png" % (save_directory, version, frame)
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

