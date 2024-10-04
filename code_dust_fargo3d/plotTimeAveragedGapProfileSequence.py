"""
plot azimuthally averaged density
then, makes movies

Usage:
python plotAveragedDensity.py frame_number <== plot one frame
python plotAveragedDensity.py -m <== make a movie instead of plots (plots already exist)
python plotAveragedDensity.py -1 <<<===== Plots a sample
python plotAveragedDensity.py <== plot all frames and make a movie

"""

import sys, os, subprocess, csv
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
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
    parser.add_argument('--dir', dest = "save_directory", default = "timeAveragedGapProfileSequences",
                         help = 'save directory (default: timeAveragedGapProfile)')
    parser.add_argument('--mpi', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpiio output files (default: do not use mpi)')
    parser.add_argument('--merge', dest = "merge", type = int, default = 0,
                         help = 'number of cores needed to merge data outputs (default: 0)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'maximum density (default: 1.1 times the max)')
    parser.add_argument('--y2_range', dest = "y2_range", type = float, nargs = 2, default = [0, 0.025],
                         help = 'range in y-axis (default: [-0.2, 0.2])')

    parser.add_argument('-l', dest = "maximum_condition", action = 'store_true', default = False,
                         help = 'plot maximum condition (default: do not do it!)')

    parser.add_argument('--derivative', dest = "derivative", action = 'store_true', default = False,
                         help = 'show derivative (default: do not do it!)')

    parser.add_argument('--zero', dest = "zero", action = 'store_true', default = False,
                         help = 'plot density at t = 0 for reference (default: do not do it!)')

    parser.add_argument('--log', dest = "log", action = 'store_true', default = False,
                         help = 'plot log y-axis (default: do not do it!)')

    parser.add_argument('-r', dest = "ref", type = int, default = 0,
                         help = 'ref case to compare to (num Jup) (default: 0)')

    parser.add_argument('--hill', dest = "hill", type = int, default = 1,
                         help = 'number of hill radii to exclude (default: 1)')

    parser.add_argument('--compare', dest = "compare", nargs = '+', default = None,
                         help = 'compare to another directory (default: do not do it!)')

    # Quantity to plot for maximum condition
    parser.add_argument('--rossby', dest = "rossby", action = 'store_true', default = False,
                         help = 'plot rossby number instead of vorticity (default: plot vorticity)')
    parser.add_argument('--residual', dest = "residual", action = 'store_true', default = False,
                         help = 'use v_theta or v_theta - v_kep (default: do not use residual)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 21,
                         help = 'fontsize of plot annotations (default: 21)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 4,
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
surface_density_power = p.sigmaslope

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
frame_range = args.frames

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
max_y = args.max_y
y2_range = args.y2_range

# Quantity to Plot
rossby = args.rossby
residual = args.residual

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

# Planet File
# Data
data = np.loadtxt("planet0.dat")
times = data[:, 0]; base_mass = data[:, 7]
accreted_mass = data[:, 8] / jupiter_mass
omegas = data[:, 10]

planet_x = data[:, 1]
planet_y = data[:, 2]

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

##### HELPER FUNCTION #####

def get_rid_of_hill(density, angle = np.pi, num_hill = 1.0):
   planet_r = 1.0
   planet_phi = angle

   x_p = planet_r * np.cos(planet_phi)
   y_p = planet_r * np.sin(planet_phi)

   hill_r = planet_r * np.power(planet_mass * jupiter_mass / 3.0, 1.0 / 3.0)
   start_r = planet_r - num_hill * hill_r
   end_r = planet_r + num_hill * hill_r

   start_ri = np.searchsorted(rad, start_r) - 1
   end_ri = np.searchsorted(rad, end_r) + 1

   for r_i, r in enumerate(rad[start_ri : end_ri]):
       x_i = r * np.cos(theta)
       y_i = r * np.sin(theta)

       distance = np.sqrt(np.power(x_i - x_p, 2) + np.power(y_i - y_p, 2))
       mask = np.ones(len(theta))
       mask[distance < num_hill * hill_r] = 0 # positive is safe, negative should be masked

       density[start_ri + r_i] *= mask

   return density

###############################################################################

##### PLOTTING #####

alpha = 0.7

labelsize = 18
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    for i, frame in enumerate(frame_range):
        density = fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)
        density_plus_one = fromfile("gasdens%d.dat" % (frame + 1)).reshape(num_rad, num_theta)
        density = (density + density_plus_one) / 2.0

        if args.hill > 0:
          density = get_rid_of_hill(density, num_hill = args.hill)
        averagedDensity = np.average(density, axis = 1)

        density_zero = fromfile("gasdens0.dat").reshape(num_rad, num_theta)
        averagedDensity_zero = np.average(density_zero, axis = 1)

        normalized_density = averagedDensity / averagedDensity_zero

        ### Plot ###
        x = rad
        y = normalized_density
        result,  = plot.plot(x, y, linewidth = linewidth, c = "b", zorder = 99, label = r"$t$ $=$ $%d$ $T_\mathrm{p}$" % frame)

    # Reference
    if args.ref > 0:
      colors = ['k', 'grey']
      #times = [20, 30, 50, 100, 140]
      times = [140]
      for i, time_i in enumerate(times):
         ref_density = []
         with open("Mt%dAm3-t%d.csv" % (args.ref, time_i), "r") as f:
            reader = csv.reader(f)
            for row in reader:
               ref_density.append(row)
         ref_density = np.array(ref_density).astype(np.float)

         #print np.shape(ref_density), ref_density

         xref = ref_density[:, 0] / 6.0 # Rp = 6 in Aoyama+Bai 23
         yref = ref_density[:, 1]
         plot.plot(xref, yref, linewidth = linewidth, c = colors[i % 2], zorder = 2, label = "A&B 2023 (t = 140)")

    plot.legend(loc = "lower right", fontsize = fontsize - 5)

    this_x = planet_x[frame]
    this_y = planet_y[frame]
    planet_location = np.sqrt(np.power(this_x, 2) + np.power(this_y, 2))

    if not args.log:
        plot.scatter([planet_location], [1.05 * 10**(-3)], c = 'k', s = 75, alpha = 0.8, clip_on = False)
    else:
        plot.scatter([planet_location], [1.005 * 10**(-2)], c = 'k', s = 75, alpha = 0.8, clip_on = False)

    if args.zero:
        #density_zero = fromfile("gasdens0.dat").reshape(num_rad, num_theta)
        #averagedDensity_zero = np.average(density_zero, axis = 1)
        normalized_density_zero = averagedDensity_zero / surface_density_zero

        x = rad
        y_zero = normalized_density_zero
        result = plot.plot(x, y_zero, linewidth = linewidth, zorder = 0)

    if args.compare is not None:
        directories = args.compare
        for i, directory in enumerate(directories):
            density_compare = (fromfile("%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta))
            density_compare_plus_one = (fromfile("%s/gasdens%d.dat" % (directory, frame + 1)).reshape(num_rad, num_theta))
            density_compare = (density_compare + density_compare_plus_one) / 2.0
            averagedDensity_compare = np.average(density_compare, axis = 1)

            if args.hill > 0:
               density_compare = get_rid_of_hill(density_compare, num_hill = args.hill)

            density_compare_zero = (fromfile("%s/gasdens0.dat" % (directory)).reshape(num_rad, num_theta))
            averagedDensity_compare_zero = np.average(density_compare_zero, axis = 1)

            normalized_density_compare = averagedDensity_compare / averagedDensity_compare_zero

            ### Plot ###
            x = rad
            y_compare = normalized_density_compare
            result = plot.plot(x, y_compare, linewidth = linewidth, alpha = 0.6, zorder = 99, label = directory)

        plot.legend(loc = "upper right")

    if args.derivative:
        twin = ax.twinx()

        ### Plot ###
        dr = rad[1] - rad[0]
        normalized_density_derivative = np.diff(normalized_density) / dr

        x2 = rad[1:]
        y2 = normalized_density_derivative
        result = twin.plot(x2, y2, c = "purple", linewidth = linewidth, alpha = 0.6, zorder = 99, label = "derivative")

        plot.legend()

        twin.set_ylim(-10, 10)

    # Axes
    if args.max_y is None:
        x_min_i = np.searchsorted(x, x_min)
        x_max_i = np.searchsorted(x, x_max)
        max_y = 1.1 * max(y[x_min_i : x_max_i])
    else:
        max_y = args.max_y

    ax.set_xlim(x_min, x_max)

    if args.log:
        ax.set_yscale("log")
        ax.set_ylim(1.0e-2, max_y)
    else:
        ax.set_ylim(0, max_y)

    # Annotate Axes
    orbit = (dt / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    current_mass += accreted_mass[frame]

    #title = readTitle()

    unit = "r_\mathrm{p}"
    ax.set_xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
    ax.set_ylabel(r"$\Sigma$ $/$ $\Sigma_0$", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    if args.log:
        y_text = 2.15

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)
    title1 = r"$h/r = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    title2 = r"$M_\mathrm{p}\ =\ %.2f$ $M_\mathrm{Jup}$" % (current_mass)
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)
    #ax.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')

    if args.maximum_condition:
        bump, _ = az.get_radial_peak(normalized_density, fargo_par, end = 1.6)
        plot.plot([bump, bump], [0, max_y], c = "b", linewidth = linewidth - 2, alpha = alpha, linestyle = "--", zorder = 20)

        twin = ax.twinx()

        vrad = (fromfile("gasvy%d.dat" % frame).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
        vtheta = (fromfile("gasvx%d.dat" % frame).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
        vorticity = utilVorticity.velocity_curl(vrad, vtheta, rad, theta, rossby = rossby, residual = residual, omega = omegas[frame])

        averaged_vorticity = np.average(vorticity, axis = 1)
        #averaged_density = np.average(normalized_density, axis = 1) # normalized_density
        maximum_condition = (normalized_density[1:] / averaged_vorticity) * (np.power(scale_height, 2) / np.power(rad[1:], 1))

        x2 = rad[1:]
        y2 = maximum_condition
        result2, = twin.plot(x2, y2, c = 'darkviolet', linewidth = linewidth, zorder = 99)

        bump, _ = az.get_radial_peak(maximum_condition, fargo_par, end = 1.6)
        plot.plot([bump, bump], y2_range, c = "darkviolet", linewidth = linewidth - 2, alpha = alpha, linestyle = "--", zorder = 20)

        # Axes
        twin.set_ylim(y2_range[0], y2_range[1])
        twin.set_yticks(np.arange(y2_range[0], y2_range[1] + 1e-9, 0.005))

        twin.set_ylabel(r"$\Sigma$ $/$ ($\nabla \times v$)$_\mathrm{z}$", fontsize = fontsize, rotation = 270, labelpad = 25)

        tkw = dict(size=4, width=1.5)
        ax.tick_params(axis = 'y', colors = result.get_color(), **tkw)
        twin.tick_params(axis = 'y', colors = result2.get_color(), **tkw)

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1].split("-")[0]
    frame_string = ""
    for frame in frame_range:
        frame_string += ("%04d-" % frame)
    frame_string = frame_string[:-1] # get rid of trailing '-'

    if version is None:
        save_fn = "%s/timeAveragedGapProfile_%s_%s.png" % (save_directory, directory_name, frame_string)
    else:
        save_fn = "%s/v%04d_timeAveragedGapProfile_%s_%s.png" % (save_directory, version, directory_name, frame_string)
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

