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
from scipy.optimize import curve_fit

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

    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "contrasts",
                         help = 'save directory (default: .)')

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

    parser.add_argument('--min_mass', dest = "min_mass", type = float, default = 0.1,
                         help = 'minimum mass on plot (default: 0.1 Jupiter mass)')
    parser.add_argument('--max_mass', dest = "max_mass", type = float, default = 1.0,
                         help = 'maximum mass on plot (default: 1.0 Jupiter mass)')
    parser.add_argument('--delta_mass', dest = "delta_mass", type = float, default = 0.1,
                         help = 'delta mass on plot (default: 0.1 Jupiter mass)')


    parser.add_argument('--negative', dest = "negative", action = 'store_true', default = False,
                         help = 'add negative mass (default: do not)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 17,
                         help = 'fontsize of plot annotations (default: 17)')
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
    shifted_density = np.roll(normalized_density, shift_c, axis = -1)
    return shifted_density, shift_c

###############################################################################

### Data ###

def get_contrasts(args_here):
    # Unwrap Args
    i, frame, directory, comparison_directory = args_here

    # Get Data
    density = fromfile("%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta) / surface_density_zero
    density, shift_c = shift_density(density, fargo_par, reference_density = density)
    azimuthal_profile = az.get_mean_azimuthal_profile(density, fargo_par, sliver_width = 1.5)

    maxima_over_time[i] = np.percentile(azimuthal_profile, 90)
    minima_over_time[i] = np.percentile(azimuthal_profile, 5)
    differences_over_time[i] = maxima_over_time[i] - minima_over_time[i]
    contrasts_over_time[i] = maxima_over_time[i] / minima_over_time[i]

    print i, frame, contrasts_over_time[i], maxima_over_time[i], minima_over_time[i], differences_over_time[i]

    # Do Comparison
    if frame >= 2500 and frame < 3000:
        density = fromfile("%s/gasdens%d.dat" % (comparison_directory, frame)).reshape(num_rad, num_theta) / surface_density_zero
        density, shift_c = shift_density(density, fargo_par, reference_density = density)
        azimuthal_profile = az.get_mean_azimuthal_profile(density, fargo_par, sliver_width = 1.5)

        maxima = np.percentile(azimuthal_profile, 90)
        minima = np.percentile(azimuthal_profile, 5)
        contrasts_over_time_comparison[i] = maxima / minima

        print i, frame, contrasts_over_time_comparison[i]


## Use These Frames ##
rate = 20 # 5 works better, but is very slow
start = 2500
max_frame = 3000 #util.find_max_frame()
frame_range_comparison = np.array(range(start, max_frame + 1, rate))

maxima_over_time = mp_array("d", len(frame_range))
minima_over_time = mp_array("d", len(frame_range))
differences_over_time = mp_array("d", len(frame_range))
contrasts_over_time = mp_array("d", len(frame_range))

contrasts_over_time_comparison = mp_array("d", len(frame_range))

for i, frame in enumerate(frame_range):
    get_contrasts((i, frame, ".", "../hi_res_high_density-2000-m2500"))

pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

#p = Pool(num_cores)
#p.map(get_extents, pool_args)
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
    # Figure
    fig = plot.figure(figsize = (7, 4), dpi = dpi)
    ax = fig.add_subplot(111)

    # Plot
    x = frame_range
    y = contrasts_over_time

    plot.plot(x, y, c = 'k', linewidth = linewidth, label = "Default")

    # Fit Data
    def exponential_decay(x, a, b):
        return 1.0 + a * np.exp(-(x - 2500.0) / b)

    comp_start = np.searchsorted(frame_range, 2500)
    comp_end = np.searchsorted(frame_range, 3000)

    popt, pcov = curve_fit(exponential_decay, x[comp_start:], y[comp_start:], p0 = [0.5, 2000])
    x_fit = np.array(range(2500, 8001, 1))
    #y_fit = fit[0] * np.power(x_fit, 2) + fit[1] * x_fit + fit[2]
    y_fit = exponential_decay(x_fit, popt[0], popt[1])

    fit = plot.plot(x_fit, y_fit, c = 'r', linewidth = linewidth, linestyle = "--", label = "Fit (default)")

    # Plot comparison
    x_comp = frame_range[comp_start : comp_end]
    y_comp = contrasts_over_time_comparison[comp_start : comp_end]

    plot.plot(x_comp, y_comp, c = 'cornflowerblue', linewidth = linewidth - 1, label = r"Restart w/ $M_\mathrm{p} = 0.05$ $M_\mathrm{J}$", zorder = 90)

    popt_comp, pcov_comp = curve_fit(exponential_decay, x_comp, y_comp, p0 = [0.5, 2000])
    x_fit_compare = np.array(range(2500, 5001, 1))
    y_fit_compare = exponential_decay(x_fit_compare, popt_comp[0], popt_comp[1])

    fit = plot.plot(x_fit_compare, y_fit_compare, c = 'orange', linewidth = linewidth - 1, linestyle = "--", label = "Fit (restart)")

    # Plot Reference (Assume vortex dies at Contrast = 1.05)
    plot.plot([x[0], x_fit[-1]], [1.05, 1.05], linewidth = 1, c = "k")
    
    # Axes
    plot.xlim(x[0], x_fit[-1])
    plot.ylim(1, 1.05 * max(y))

    # Annotate
    plot.xlabel("Time (planet orbits)", fontsize = fontsize)
    #plot.ylabel(r"Contrast ($C$ $\equiv$ $\Sigma_\mathrm{peak}$ $/$ $\Sigma_\mathrm{trough}$)", fontsize = fontsize)
    plot.ylabel(r"Azimuthal Contrast", fontsize = fontsize)

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    title1 = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $M_\mathrm{p}=%.2f$ $M_\mathrm{J}$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, planet_mass)
    #plot.title("%s" % (title1), y = 1.01, fontsize = fontsize + 1)

    plot.legend(loc = "upper right", fontsize = fontsize - 6)

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1].split("-")[0]

    if version is None:
        save_fn = "%s/%s_contrastComparisonOverTime.png" % (save_directory, directory_name)
    else:
        save_fn = "%s/v%04d_%s_contrastComparisonOverTime.png" % (save_directory, version, directory_name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)