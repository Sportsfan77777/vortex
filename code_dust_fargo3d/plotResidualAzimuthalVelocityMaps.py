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
from readTitle import readTitle

from advanced import Parameters
from reader import Fields

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
    parser.add_argument('--dir', dest = "save_directory", default = "radialVelocityMaps",
                         help = 'save directory (default: radialVelocityMaps)')
    parser.add_argument('--mpi', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpio output files (default: use dat)')
    parser.add_argument('--merge', dest = "merge", type = int, default = 0,
                         help = 'number of cores needed to merge data outputs (default: 0)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--shift', dest = "center", action = 'store_true', default = False,
                         help = 'center frame on vortex peak or middle (default: do not center)')

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
    parser.add_argument('--cmap', dest = "cmap", default = "seismic",
                         help = 'color map (default: seismic)')
    parser.add_argument('--cmax', dest = "cmax", type = float, default = 0.2,
                         help = 'maximum radial velocity in colorbar (default: 0.2)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
p = Parameters()

num_rad = p.ny; num_theta = p.nx
r_min = p.ymin; r_max = p.ymax

surface_density_zero = p.sigma0

planet_mass = 1.0
taper_time = p.masstaper

dt = p.ninterm * p.dt

"""
fargo_par = util.get_pickled_parameters()

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

size = fargo_par["PSIZE"]
"""

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

merge = args.merge
mpi = args.mpi

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
center = args.center

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
clim = [-args.cmax, args.cmax]

fontsize = args.fontsize
dpi = args.dpi

"""
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
center = args.center

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
"""

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

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
    if merge > 0:
        num_merged_cores = merge
        gas_density = util.read_merged_data(frame, num_merged_cores, num_rad, num_theta)
        velocity = util.read_merged_data(frame, num_merged_cores, num_rad, num_theta, field = "vx")
    elif mpi:
        gas_density = Fields("./", 'gas', frame).get_field("dens").reshape(num_rad, num_theta)
        velocity = Fields("./", 'gas', frame).get_field("vx").reshape(num_rad, num_theta)
    else:
        gas_density = fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)
        velocity = fromfile("gasvx%d.dat" % frame).reshape(num_rad, num_theta)
    normalized_gas_density = gas_density / surface_density_zero

    # Take Residual
    keplerian_velocity = rad * (np.power(rad, -1.5) - 1) # in rotating frame, v_k = r * (r^-1.5 - r_p^-1.5)
    azimuthal_velocity -= keplerian_velocity[:, None]

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi)
    result = ax.pcolormesh(x, y, np.transpose(velocity), cmap = cmap)

    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    if use_contours:
        levels = np.linspace(low_contour, high_contour, num_levels)
        colors = generate_colors(num_levels)
        plot.contour(x, y, np.transpose(normalized_density), levels = levels, origin = 'upper', linewidths = 1, colors = colors)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)

    # Annotate Axes
    orbit = (dt / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    #title = readTitle()

    unit = "r_\mathrm{p}"
    plot.xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)

    title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/radialVelocityMap_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_radialVelocityMap_%04d.png" % (save_directory, version, frame)
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

