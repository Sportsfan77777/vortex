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
from matplotlib.ticker import MultipleLocator

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

problem_frames = np.array([3000, 3125, 3200, 3250, 3550, 7425]) # h = 0.08, A = 0.167

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "rossbyWavesOverTime",
                         help = 'save directory (default: .)')

    # Quantity to plot
    parser.add_argument('--rossby', dest = "rossby", action = 'store_true', default = False,
                         help = 'plot rossby number instead of vorticity (default: plot vorticity)')
    parser.add_argument('--residual', dest = "residual", action = 'store_true', default = False,
                         help = 'use v_theta or v_theta - v_kep (default: do not use residual)')

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

    parser.add_argument('--min', dest = "min_mass", type = float, default = 0.1,
                         help = 'minimum mass on plot (default: 0.1 Jupiter mass)')
    parser.add_argument('--max', dest = "max_mass", type = float, default = 1.0,
                         help = 'maximum mass on plot (default: 1.0 Jupiter mass)')
    parser.add_argument('--delta', dest = "delta_mass", type = float, default = 0.1,
                         help = 'delta mass on plot (default: 0.1 Jupiter mass)')
    parser.add_argument('--minor_delta', dest = "minor_delta_mass", type = float, default = None,
                         help = 'delta mass on plot (default: 0.1 Jupiter mass)')

    parser.add_argument('-r', dest = "check_rossby", type = int, default = 1000000,
                         help = 'frame at which you start using the Rossby number for measuring everything (default: infinity)')
    parser.add_argument('-e', dest = "extreme_cutoff", type = int, default = 1000000,
                         help = 'frame at which you start using the extreme Rossby number cutoff (default: infinity)')
    parser.add_argument('-i', dest = "include_aspect", action = 'store_true', default = False,
                         help = 'include aspect ratio (default: do not)')
    parser.add_argument('--negative', dest = "negative", action = 'store_true', default = False,
                         help = 'add negative mass (default: do not)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 17,
                         help = 'fontsize of plot annotations (default: 17)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 15,
                         help = 'fontsize of plot annotations (default: 15)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 2,
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

scale_height = p.aspectratio
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

# Quantity to Plot
rossby = args.rossby
residual = args.residual

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

check_rossby = args.check_rossby
extreme_cutoff = args.extreme_cutoff
negative = args.negative

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

# Planet Data
data = np.loadtxt("planet0.dat")
times = data[:, 0]
base_mass = data[:, 7]
accreted_mass = data[:, 8]

total_mass = base_mass + accreted_mass

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Helper Functions ###

def shift_density(normalized_density, vorticity, fargo_par, option = "away", reference_density = None, frame = None):
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
    shifted_vorticity = np.roll(vorticity, shift_c)
    shifted_density = np.roll(normalized_density, shift_c, axis = -1)
    return shifted_density, shifted_vorticity, shift_c

###############################################################################

### Data ###

def get_rossby_criteria(args_here):
    # Unwrap Args
    i, frame = args_here

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

    # Diagnostics
    peak_rad, peak_density = az.get_radial_peak(averaged_density, fargo_par)
    peak_rad_i = np.searchsorted(rad, peak_rad)
    inner_limit_i = np.searchsorted(rad, 1.1) # Make inner limit a variable in the future
    outer_limit_i = np.searchsorted(rad, 2.0) # Make outer limit a variable in the future

    inner_max_diff_i = np.argmax(diff_maximum_condition[inner_limit_i : peak_rad_i])
    outer_max_diff_i = np.argmin(diff_maximum_condition[peak_rad_i : outer_limit_i])

    inner_max_diff_i += inner_limit_i # put the "inner disk" back
    outer_max_diff_i += peak_rad_i

    inner_rossby_rad = rad[inner_max_diff_i]
    outer_rossby_rad = rad[outer_max_diff_i]
    difference = outer_rossby_rad - inner_rossby_rad

    inner_rossby_value = diff_maximum_condition[inner_max_diff_i]
    outer_rossby_value = diff_maximum_condition[outer_max_diff_i] * -1.0 # absolute value

    # Store Data
    inner_rossby_rad_over_time[i] = inner_rossby_rad
    peak_rad_over_time[i] = peak_rad
    outer_rossby_rad_over_time[i] = outer_rossby_rad

    rossby_rad_difference_over_time[i] = outer_rossby_rad - inner_rossby_rad
    inner_peak_difference_over_time[i] = peak_rad - inner_rossby_rad
    outer_peak_difference_over_time[i] = outer_rossby_rad - peak_rad

    inner_rossby_value_over_time[i] = inner_rossby_value
    outer_rossby_value_over_time[i] = outer_rossby_value

    print i, frame, "Rad: ", inner_rossby_rad_over_time[i], peak_rad_over_time[i], outer_rossby_rad_over_time[i]
    print i, frame, "Diff:", rossby_rad_difference_over_time[i], inner_peak_difference_over_time[i], outer_peak_difference_over_time[i]
    print i, frame, "Val: ", inner_rossby_value_over_time[i], outer_rossby_value_over_time[i]
    #print i, frame, azimuthal_extent_over_time[i], radial_extent_over_time[i], radial_peak_over_time[i], radial_peak_over_time_a[i], contrasts_over_time[i]


## Use These Frames ##
rate = 1 # 5 works better, but is very slow
start = 50
max_frame = 100 #util.find_max_frame()
#frame_range = np.array(range(start, max_frame + 1, rate))

inner_rossby_rad_over_time = mp_array("d", len(frame_range))
outer_rossby_rad_over_time = mp_array("d", len(frame_range))
peak_rad_over_time = mp_array("d", len(frame_range))

rossby_rad_difference_over_time = mp_array("d", len(frame_range))
inner_peak_difference_over_time = mp_array("d", len(frame_range))
outer_peak_difference_over_time = mp_array("d", len(frame_range))

inner_rossby_value_over_time = mp_array("d", len(frame_range))
outer_rossby_value_over_time = mp_array("d", len(frame_range))

for i, frame in enumerate(frame_range):
    get_rossby_criteria((i, frame))

pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

#p = Pool(num_cores)
#p.map(get_rossby_criteria, pool_args)
#p.terminate()

##### Helper Functions #####

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

###############################################################################

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (8, 10), dpi = dpi)

    ##### PLOT 1 #####
    par1 = plot.subplot(3, 1, 1)

    if args.include_aspect:
        par4 = host.twinx()

    # Plot
    x = frame_range
    y1 = inner_rossby_value_over_time
    y2 = outer_rossby_value_over_time

    p1, = par1.plot(x, y1, c = 'b', linewidth = linewidth)
    p2, = par1.plot(x, y2, c = 'orange', linewidth = linewidth)

    # Axes
    #host.xaxis.set_minor_locator(MultipleLocator(500))
    par1.set_xlim(x[0], x[-1])
    par1.set_ylim(10**(-2), 10**2)

    par1.set_yscale('log')

    #par1.set_xlabel("Time (planet orbits)", fontsize = fontsize)
    par1.set_ylabel("Derivatives", fontsize = fontsize)

    ##### PLOT 2 #####
    par2 = plot.subplot(3, 1, 2)

    y1 = outer_rossby_rad_over_time
    y2 = peak_rad_over_time
    y3 = inner_rossby_rad_over_time

    p1, = par2.plot(x, y1, c = 'b', linewidth = linewidth)
    p2, = par2.plot(x, y2, c = 'orange', linewidth = linewidth)
    p3, = par2.plot(x, y3, c = 'g', linewidth = linewidth)

    par2.set_xlim(x[0], x[-1])
    par2.set_ylim(1, 2)

    par2.set_ylabel("Radii", fontsize = fontsize)

    ##### PLOT 3 #####
    par3 = plot.subplot(3, 1, 3)

    y1 = rossby_rad_difference_over_time
    y2 = inner_peak_difference_over_time
    y3 = outer_peak_difference_over_time

    p1, = par3.plot(x, y1, c = 'b', linewidth = linewidth)
    p2, = par3.plot(x, y2, c = 'orange', linewidth = linewidth)
    p3, = par3.plot(x, y3, c = 'g', linewidth = linewidth)

    par3.set_xlim(x[0], x[-1])
    par3.set_ylim(0, 0.5)

    par3.set_xlabel("Time (planet orbits)", fontsize = fontsize)
    par3.set_ylabel("Differences", fontsize = fontsize)

    # Annoatate

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    title1 = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    plot.title("%s" % (title1), y = 1.035, fontsize = fontsize + 1)

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1].split("-")[0]

    if version is None:
        save_fn = "%s/%s_rossbyWavesOverTime.png" % (save_directory, directory_name)
    else:
        save_fn = "%s/v%04d_%s_rossbyWavesOverTime.png" % (save_directory, version, directory_name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)