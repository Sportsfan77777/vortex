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
    parser.add_argument('--dir', dest = "save_directory", default = "dustVariance",
                         help = 'save directory (default: gasDensityMaps)')
    parser.add_argument('--mpi', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpio output files (default: use dat)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'maximum density (default: 1.1 times the max)')

    parser.add_argument('--zero', dest = "zero", action = 'store_true', default = False,
                         help = 'plot density at t = 0 for reference (default: do not do it!)')

    parser.add_argument('--compare', dest = "compare", default = None,
                         help = 'compare to fargo (default: do not do it!)')
    parser.add_argument('--data', dest = "data", default = None,
                         help = 'compare to data from another directory (default: do not do it!)')
    
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
dust_surface_density_zero = p.sigma0 * p.epsilon
taper_time = p.masstaper

viscosity = p.nu
scale_height = p.aspectratio

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

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
#fargo_par["rad"] = rad
#fargo_par["theta"] = theta

###############################################################################

### Helper Functions ###

def find_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max radial density before r = 2.3
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start : outer_disk_end])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_density

### Data ###

def get_excess_mass(args_here):
    # Unwrap Args
    i, frame = args_here

    # Get Data
    if mpi:
        field = "dens"
        density = Fields("./", 'gas', frame).get_field(field).reshape(num_rad, num_theta) / surface_density_zero
        #background_density = Fields("./", 'gas', frame - 1).get_field(field).reshape(num_rad, num_theta) / surface_density_zero
    else:
        density = fromfile("dust1dens%d.dat" % frame).reshape(num_rad, num_theta) / dust_surface_density_zero
        #background_density = fromfile("dust1dens%d.dat" % (frame - 1)).reshape(num_rad, num_theta) / dust_surface_density_zero

    if args.compare:
        fargo_directory = args.compare
        density_compare = (fromfile("%s/dust1dens%d.dat" % (fargo_directory, frame)).reshape(num_rad, num_theta)) / dust_surface_density_zero
        #background_density_compare = (fromfile("%s/dust1dens%d.dat" % (fargo_directory, frame - 1)).reshape(num_rad, num_theta)) / dust_surface_density_zero

    def helper(density):
        diff_density = density # - background_density
        #diff_density[diff_density < 0] = 0 # only include excess

        # Extract Near Vortex
        averagedDensity = np.average(density, axis = 1)
        peak_rad, peak_density = find_peak(averagedDensity)

        vortex_start = np.max([1.0, peak_rad - 5.0 * scale_height])
        vortex_end = peak_rad + 5.0 * scale_height

        vortex_start_i = np.searchsorted(rad, vortex_start)
        vortex_end_i = np.searchsorted(rad, vortex_end)

        vortex_rad = rad[vortex_start_i : vortex_end_i]
        vortex_diff_density = diff_density[vortex_start_i : vortex_end_i]

        vortex_excess = np.average(vortex_diff_density, axis = 0)

        variance = np.std(vortex_excess)
        average = np.average(vortex_excess)

        # Add up mass
        #dr = rad[1] - rad[0] # assumes arithmetic grid
        #d_phi = theta[1] - theta[0]

        #excess_mass = np.sum((dr * d_phi) * vortex_rad[:, None] * vortex_diff_density)
        #return excess_mass, vortex_excess

        return variance / average, vortex_excess

    excess_mass, vortex_excess = helper(density)
    if args.compare:
        excess_mass_compare, vortex_excess_compare = helper(density_compare)

    # Get Peak
    peak_diff_density = np.max(vortex_excess)
    if args.compare:
        peak_diff_density_compare = np.max(vortex_excess_compare)

    # Print Update
    print "%d: %.4f, %.4f" % (frame, excess_mass, peak_diff_density)
    if args.compare:
        print "%d: %.4f, %.4f" % (frame, excess_mass_compare, peak_diff_density_compare)

    # Store Data
    mass_over_time[i] = excess_mass
    peak_over_time[i] = peak_diff_density

    if args.compare:
        mass_over_time_compare[i] = excess_mass_compare
        peak_over_time_compare[i] = peak_diff_density_compare

###############################################################################

## Use These Frames ##
rate = 1 # 5 works better, but is very slow
start = 50
max_frame = 100 #util.find_max_frame()
#frame_range = np.array(range(start, max_frame + 1, rate))

#mass_over_time = np.zeros(len(frame_range))
#peak_over_time = np.zeros(len(frame_range))

mass_over_time = mp_array("d", len(frame_range))
peak_over_time = mp_array("d", len(frame_range))

mass_over_time_compare = mp_array("d", len(frame_range))
peak_over_time_compare = mp_array("d", len(frame_range))

#for i, frame in enumerate(frame_range):
#    get_excess_mass((i, frame))

pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

p = Pool(num_cores)
p.map(get_excess_mass, pool_args)
p.terminate()

max_mass = np.max(mass_over_time)
max_peak = np.max(peak_over_time)

if args.compare:
    max_mass_compare = np.max(mass_over_time_compare)
    max_peak_compare = np.max(peak_over_time_compare)

## Pickle to combine later ##

pickle.dump(np.array(frame_range), open("dust_variance_frames.p", "wb"))
pickle.dump(np.array(mass_over_time), open("dust_variance_values.p", "wb"))

##### PLOTTING #####

def make_plot(show = False):
    # Figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)

    # Curves
    plot.plot(frame_range, mass_over_time, linewidth = linewidth)
    #plot.plot(frame_range, peak_over_time, linewidth = linewidth - 1, label = "Peak")
    if args.compare:
        plot.plot(frame_range, mass_over_time_compare, linewidth = linewidth, label = "compare")

    if args.data:
        frame_range_data = pickle.load(open("%s/dust_variance_frames.p" % args.data, "rb"))
        mass_over_time_data = pickle.load(open("%s/dust_variance_frames.p" % args.data, "rb"))
        plot.plot(frame_range_data, mass_over_time_data, linewidth = linewidth, label = "data")

    # Reference Lines
    plot.plot([0, frame_range[-1]], 0.10 * np.ones(2), linewidth = 2, color = "black")
    #plot.plot([0, frame_range[-1]], 0.10 * max_mass * np.ones(2), linewidth = 2, color = "black")
    #plot.plot([0, frame_range[-1]], 0.10 * max_peak * np.ones(2), linewidth = 1, color = "black")
    if args.compare:
        plot.plot([0, frame_range[-1]], 0.10 * max_mass_compare * np.ones(2), linewidth = 2, color = "black")

    # Annotate
    #this_title = readTitle()
    title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Dust Variation", fontsize = fontsize)

    title1 = os.getcwd().split("/")[-1]
    plot.title(title1, fontsize = fontsize)

    #plot.legend(loc = "upper left")

    # Limits
    plot.xlim(frame_range[0], frame_range[-1])
    #plot.ylim(0.0, 1.0)

    plot.yscale('log')

    # Save + Close
    directory_name = os.getcwd().split("/")[-1].split("-")[0]

    if version is None:
        save_fn = "%s/%s_dustVariation.png" % (save_directory, directory_name)
    else:
        save_fn = "%s/v%04d_%s_dustVariation.png" % (save_directory, version, directory_name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
