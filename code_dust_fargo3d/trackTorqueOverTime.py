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
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np
from scipy import signal as sig
from scipy.ndimage import filters as ff

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
import azimuthal as az
#from readTitle import readTitle

from advanced import Parameters
#from reader import Fields

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
    parser.add_argument('--dir', dest = "save_directory", default = "torque",
                         help = 'save directory (default: torque)')

    # Reference
    parser.add_argument('--ref', dest = "ref", type = int, default = 0,
                         help = 'reference taper time for prescribed growth curve (default: no reference)')
    parser.add_argument('--compare', dest = "compare", nargs = '+', default = None,
                         help = 'select directories to compare planet growth rates')


    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'maximum density (default: 1.1 times the max)')

    parser.add_argument('--offset', dest = "offset", type = float, default = 0.0,
                         help = 'time offset for compare (default: 0.0)')

    parser.add_argument('--negative', dest = "negative", action = 'store_true', default = False,
                         help = 'add negative mass (default: do not)')
    
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

taper_time = p.masstaper

viscosity = p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]


"""
num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]


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

# Reference
ref = args.ref

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
if args.r_lim is None:
    x_min = 0; x_max = 1000
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
max_y = args.max_y

negative = args.negative

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

# Planet
data = np.loadtxt("planet0.dat")
times = data[:, 0]
planet_x = data[:, 1]
planet_y = data[:, 2]
planet_radii = np.sqrt(np.power(planet_x, 2) + np.power(planet_y, 2))

BigG = 1.0

### Add new parameters to dictionary ###
#fargo_par["rad"] = rad
#fargo_par["theta"] = theta

### Helper Functions ###

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
ks = 40.0 # Kernel Size
ks_small = ks / 5.0 # Smaller kernel to check the normal kernel

def get_torque(args_here):
    # Unwrap Args
    i, frame = args_here

    # Get Data
    density = fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)

    if args.compare:
        fargo_directory = args.compare
        density_compare = (fromfile("%s/gasdens%d.dat" % (fargo_directory, frame)).reshape(num_rad, num_theta))

    def helper(density):
        frame_i = np.searchsorted(times, frame)
        px = planet_x[frame_i]
        py = planet_y[frame_i]
        planet_r = planet_radii[frame_i]

        # Torque
        r_element = np.array([np.outer(rad, np.cos(theta)), np.outer(rad, np.sin(theta))]) # star to fluid element 
        r_diff = np.array([np.outer(rad, np.cos(theta)) - px, np.outer(rad, np.sin(theta)) - py]) # planet to fluid element
        dist_sq = np.einsum("ijk,ijk->jk", r_diff, r_diff)

        coeff = BigG * planet_mass * density / dist_sq
        direction = np.cross(r_element, r_diff, axis = 0)

        torque_density_per_area = coeff * direction
        d_rad = rad[1] - rad[0]; d_theta = theta[1] - theta[0]
        area = rad[:, None] * np.outer(d_rad, d_theta)

        torque_density = torque_density_per_area * area
        normalized_torque_density = torque_density / surface_density_zero # / np.sqrt(2.0 * np.pi) / scale_height_function[:, None]

        radial_torque_density_profile = np.sum(torque_density, axis = -1)

        # Split the disc
        planet_ri = np.searchsorted(rad, planet_r)
        inner_torque = np.sum(radial_torque_density_profile[:planet_ri])
        outer_torque = np.sum(radial_torque_density_profile[planet_ri:])

        return inner_torque, outer_torque

    inner_torque, outer_torque = helper(density)

    if args.compare:
        inner_torque_compare, outer_torque_compare = helper(density_compare)

    # Print Update
    print "%d: %.10f, %.10f" % (frame, inner_torque, outer_torque)
    if args.compare:
        print "%d: %.10f, %.10f" % (frame, inner_torque_compare, outer_torque_compare)

    # Store Data
    inner_torque_over_time[i] = inner_torque
    outer_torque_over_time[i] = outer_torque

    if args.compare:
        inner_torque_over_time_compare[i] = inner_torque_compare
        outer_torque_over_time_compare[i] = outer_torque_compare

###############################################################################

## Use These Frames ##
rate = 1 # 5 works better, but is very slow
start = 50
max_frame = 100 #util.find_max_frame()
#frame_range = np.array(range(start, max_frame + 1, rate))

#torque_over_time = np.zeros(len(frame_range))
#peak_over_time = np.zeros(len(frame_range))

inner_torque_over_time = mp_array("d", len(frame_range))
outer_torque_over_time = mp_array("d", len(frame_range))

if args.compare:
    inner_torque_over_time_compare = mp_array("d", len(frame_range))
    outer_torque_over_time_compare = mp_array("d", len(frame_range))

#for i, frame in enumerate(frame_range):
#    get_excess_mass((i, frame))

pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

p = Pool(num_cores)
p.map(get_torque, pool_args)
p.terminate()

if args.compare:
    torque_compare = np.max(torque_over_time_compare)

## Pickle to combine later ##

pickle.dump(np.array(frame_range), open("torque_frames.p", "wb"))
pickle.dump(np.array(inner_torque_over_time), open("inner_torque_values.p", "wb"))
pickle.dump(np.array(outer_torque_over_time), open("outer_torque_values.p", "wb"))

###############################################################################

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Data ###
    # Planet
    cwd = os.getcwd().split("/")[-1]

    ### Plot ###
    kernel_size = 5

    x = frame_range
    y1 = smooth(inner_torque_over_time, kernel_size)
    y2 = smooth(outer_torque_over_time, kernel_size)
    result1 = plot.plot(x, y1, linewidth = linewidth, zorder = 99, label = "%s (Inner)" % cwd)
    result2 = plot.plot(x, y2, linewidth = linewidth, zorder = 99, label = "%s (Outer)" % cwd)

    # Raw Data
    torque_data = np.loadtxt("tqwk0.dat")
    x_raw = torque_data[:, 0]
    y1_raw = torque_data[:, 1]
    y2_raw = torque_data[:, 2]

    y1_smooth = smooth(y1_raw, ks)
    y2_smooth = smooth(y2_raw, ks)

    raw1 = plot.plot(x_raw, y1_smooth, linewidth = linewidth, zorder = 99, label = "%s (Inner Raw)" % cwd)
    raw2 = plot.plot(x_raw, y2_smooth, linewidth = linewidth, zorder = 99, label = "%s (Outer Raw)" % cwd)

    if args.ref > 0:
        x = times
        y_ref = np.power(np.sin((np.pi / 2) * (1.0 * times / args.ref)), 2) * 1.0
        plot.plot(x, y_ref, linewidth = linewidth, alpha = 0.5)

    if args.compare is not None:
        for i, directory in enumerate(args.compare):
            data_comp = np.loadtxt("%s/planet0.dat" % directory)
            times = data_comp[:, 0] + args.offset

            planet_x = data_comp[:, 1]
            planet_y = data_comp[:, 2]
            planet_radii = np.sqrt(np.power(planet_x, 2) + np.power(planet_y, 2))

            x_comp = times
            y_comp = planet_radii
            result = plot.plot(x_comp, y_comp, linewidth = linewidth, zorder = 1, label = "%s" % directory)

    plot.legend(loc = "lower left")

    # Axes
    if args.max_y is None:
        min_y1 = -1.1 * min(y1_smooth)
        min_y2 = -1.1 * min(y2_smooth)

        max_y1 = 1.1 * max(y1_smooth)
        max_y2 = 1.1 * max(y2_smooth)

        max_y = max([min_y1, min_y2, max_y1, max_y2])
    else:
        max_y = args.max_y

    plot.xlim(x_min, x_max)
    plot.ylim(-max_y, max_y)

    #title = readTitle()

    unit = "orbits"
    plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel(r"Torque", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)

    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title1), y = 1.015, fontsize = fontsize + 1)
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')


    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1]

    if version is None:
        save_fn = "%s/%s_torqueOverTime.png" % (save_directory, directory_name)
    else:
        save_fn = "%s/v%04d_%s_torqueOverTime.png" % (save_directory, version, directory_name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
