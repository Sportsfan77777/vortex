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
    parser.add_argument('--dir', dest = "save_directory", default = "gapDepths",
                         help = 'save directory (default: gapDepths)')

    # Reference
    parser.add_argument('--ref', dest = "ref", type = int, default = 0,
                         help = 'reference taper time for prescribed growth curve (default: no reference)')
    parser.add_argument('--compare', dest = "compare", nargs = '+', default = None,
                         help = 'compare to another directory (default: do not do it!)')


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
    x_min = frame_range[0]; x_max = frame_range[-1]
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
base_mass = data[:, 7]
accreted_mass = data[:, 8] #/ jupiter_mass

BigG = 1.0

### Add new parameters to dictionary ###
#fargo_par["rad"] = rad
#fargo_par["theta"] = theta

### Helper Functions ###

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
ks = 40.0 # Kernel Size
ks_small = ks / 5.0 # Smaller kernel to check the normal kernel

def find_min(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 0.9) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 1.25) # look for max radial density before r = 2.3
    min_rad_outer_index = np.argmin(averagedDensity[outer_disk_start : outer_disk_end])

    min_index = outer_disk_start + min_rad_outer_index
    min_rad = rad[min_index]
    min_density = averagedDensity[min_index]

    return min_density

def get_min(density):
    # Get Data
    averagedDensity = np.average(density, axis = 1)
    #normalized_density = averagedDensity / surface_density_zero # FIX THIS: DIVIDE BY THE WHOLE DENSITY PROFILE!!!!!!!!!!
    
    # Get Minima
    min_density = find_min(averagedDensity)

    # Print Update

    # Store Data
    return 1.0 / min_density

def get_gap_depth(args_here):
    # Unwrap Args
    i, frame, directory = args_here

    # Get Data
    density = fromfile("%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta)
    density_zero = fromfile("%s/gasdens0.dat" % (directory)).reshape(num_rad, num_theta)
    
    gap_depth = get_min(density / density_zero)

    # Print Update
    print "%s %d: %.10f" % (directory, frame, gap_depth)

    # Store Data
    gap_depth_over_time[i] = gap_depth

###############################################################################

## Use These Frames ##
rate = 1 # 5 works better, but is very slow
start = 50
max_frame = 100 #util.find_max_frame()
#frame_range = np.array(range(start, max_frame + 1, rate))

#torque_over_time = np.zeros(len(frame_range))
#peak_over_time = np.zeros(len(frame_range))

gap_depth_over_time = mp_array("d", len(frame_range))

#for i, frame in enumerate(frame_range):
#    get_excess_mass((i, frame))

pool_args = [(i, frame, ".") for i, frame in enumerate(frame_range)]

p = Pool(num_cores)
p.map(get_gap_depth, pool_args)
p.terminate()

gap_depth_array = np.array(gap_depth_over_time)

## Pickle to combine later ##

pickle.dump(np.array(frame_range), open("gap_depth_frames.p", "wb"))
pickle.dump(np.array(gap_depth_over_time), open("gap_depth_values.p", "wb"))

###############################################################################

##### PLOTTING #####

labelsize = 18
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

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
    y = gap_depth_array
    result1 = plot.plot(x, y, c = 'b', linewidth = linewidth, zorder = 99)

    # Reference
    num_jupiters = int(round(planet_mass, 0))
    colors = ['k', 'grey']
    times = [20, 30, 50, 100, 140]
    for i, time_i in enumerate(times):
       ref_density = []
       with open("Mt%dAm3-t%d.csv" % (num_jupiters, time_i), "r") as f:
          reader = csv.reader(f)
          for row in reader:
             ref_density.append(row)
       ref_density = np.array(ref_density).astype(np.float)

       ref_density_min = 1.0 / np.min(ref_density)
       plot.scatter(time_i, ref_density_min, s = 50, c = 'k')

    # Compare
    if args.compare is not None:
        directories = args.compare
        max_yc = 0

        for i, directory in enumerate(directories):
            pool_args = [(i, frame, directory) for i, frame in enumerate(frame_range)]

            p = Pool(num_cores)
            p.map(get_gap_depth, pool_args)
            p.terminate()

            ### Plot ###
            x = frame_range
            y_compare = gap_depth_over_time
            result_compare = plot.plot(x, y_compare, linewidth = linewidth, alpha = 0.6, zorder = 99, label = directory)

            # Max
            x_min_i = np.searchsorted(x, x_min)
            x_max_i = np.searchsorted(x, x_max)
            this_max_y = 1.1 * max(y_compare[x_min_i : x_max_i])

            if this_max_y > max_yc:
                max_yc = this_max_y

        plot.legend(loc = "upper right")

    # Axes
    if args.max_y is None:
        x_min_i = np.searchsorted(x, x_min)
        x_max_i = np.searchsorted(x, x_max)
        max_y = 1.1 * max(y[x_min_i : x_max_i])

        if args.compare:
            if (max_yc > max_y):
                max_y = max_yc

    else:
        max_y = args.max_y

    plot.xlim(x_min, x_max)
    plot.ylim(0, max_y)

    #title = readTitle()

    unit = "orbits"
    plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel(r"Gap Depth", fontsize = fontsize)

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
        save_fn = "%s/%s_gapDepthOverTime.png" % (save_directory, directory_name)
    else:
        save_fn = "%s/v%04d_%s_gapDepthOverTime.png" % (save_directory, version, directory_name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
