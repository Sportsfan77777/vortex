"""
plot azimuthally averaged density
then, makes movies

Usage:
python plotAveragedDensity.py frame_number <== plot one frame
python plotAveragedDensity.py -m <== make a movie instead of plots (plots already exist)
python plotAveragedDensity.py -1 <<<===== Plots a sample
python plotAveragedDensity.py <== plot all frames and make a movie

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
import utilVorticity
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
    parser.add_argument('--dir', dest = "save_directory", default = "maximumCondition",
                         help = 'save directory (default: maximumCondition)')

    # Quantity to plot
    parser.add_argument('--rossby', dest = "rossby", action = 'store_true', default = False,
                         help = 'plot rossby number instead of vorticity (default: plot vorticity)')
    parser.add_argument('--residual', dest = "residual", action = 'store_true', default = False,
                         help = 'use v_theta or v_theta - v_kep (default: do not use residual)')
    parser.add_argument('--full', dest = "full_term", action = 'store_true', default = False,
                         help = 'use full term for extremum condition (default: just use the vortensity part)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--y_range', dest = "y_range", type = float, nargs = 2, default = [0, 0.2],
                         help = 'range in y-axis (default: [0, 0.2])')
    parser.add_argument('--y2_range', dest = "y2_range", type = float, nargs = 2, default = [-1.5, 1.5],
                         help = 'range in y-axis (default: [-1, 1])')

    parser.add_argument('--zero', dest = "zero", action = 'store_true', default = False,
                         help = 'plot density at t = 0 for reference (default: do not do it!)')

    parser.add_argument('--compare', dest = "compare", default = None,
                         help = 'compare to another directory (default: do not do it!)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'fontsize of plot annotations (default: 3)')
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

scale_height = p.aspectratio
viscosity = p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

"""
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

# Quantity to Plot
rossby = args.rossby
residual = args.residual
full_term = args.full_term

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
y_range = args.y_range
y2_range = args.y2_range

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

# Planet File
# Data
data = np.loadtxt("planet0.dat")
times = data[:, 0]; base_mass = data[:, 7]
accreted_mass = data[:, 8] / jupiter_mass

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

##### PLOTTING #####

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    host = fig.add_subplot(111)
    twin = host.twinx()

    # Data
    normalized_density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero

    vrad = (fromfile("gasvy%d.dat" % frame).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
    vtheta = (fromfile("gasvx%d.dat" % frame).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
    vorticity = utilVorticity.velocity_curl(vrad, vtheta, rad, theta, rossby = rossby, residual = residual)

    averaged_vorticity = np.average(vorticity, axis = 1)
    averaged_density = np.average(normalized_density, axis = 1)
    maximum_condition = (averaged_density[1:] / averaged_vorticity) * (np.power(scale_height, 2) / np.power(rad[1:], 1))

    dr = rad[1] - rad[0]
    diff_maximum_condition = np.diff(maximum_condition) / dr

    ### Plot ###
    x = rad[1:]
    y = maximum_condition
    y2 = diff_maximum_condition
    result = host.plot(x, y, c = 'b', linewidth = linewidth + 1, zorder = 99, label = "Condition")
    result2 = twin.plot(x[1:], y2, c = 'purple', linewidth = linewidth, zorder = 99, label = "Derivative")

    if args.zero:
        density_zero = fromfile("gasdens0.dat").reshape(num_rad, num_theta)
        averagedDensity_zero = np.average(density_zero, axis = 1)
        normalized_density_zero = averagedDensity_zero / surface_density_zero

        x = rad
        y_zero = normalized_density_zero
        result = plot.plot(x, y_zero, linewidth = linewidth, zorder = 0)

    if args.compare is not None:
        directory = args.compare
        density_compare = (fromfile("%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta))
        averagedDensity_compare = np.average(density_compare, axis = 1)
        normalized_density_compare = averagedDensity_compare / surface_density_zero

        ### Plot ###
        x = rad
        y_compare = normalized_density_compare
        result = plot.plot(x, y_compare, linewidth = linewidth, alpha = 0.6, zorder = 99, label = "compare")

    plot.legend(loc = "upper left")

    # Axes
    host.set_xlim(x_min, x_max)
    host.set_ylim(y_range[0], y_range[1])
    twin.set_ylim(y2_range[0], y2_range[1])

    # Annotate Axes
    orbit = (dt / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    current_mass += accreted_mass[frame]

    #title = readTitle()

    unit = "r_\mathrm{p}"
    host.set_xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
    host.set_ylabel(r"$c_s^2$ $\Sigma$ $/$ ($\nabla \times v$)$_\mathrm{z}$", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)
    title1 = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)
    host.text(x_mid, y_text * y_range[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')


    # Save, Show, and Close
    if version is None:
        save_fn = "%s/maximumCondition_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_maximumCondition_%04d.png" % (save_directory, version, frame)
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

