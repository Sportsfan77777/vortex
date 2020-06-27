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
from matplotlib.ticker import ScalarFormatter
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

master_directories = {}
master_directories[87] = ["h08_nu7_a167-offset", "h08_nu7_a05-offset", "h08_nu7_a02-offset", "h08_nu7_a01-offset"]
master_directories[67] = ["h06_nu7_a50-offset", "h06_nu7_a167-offset", "h06_nu7_a05-offset", "h06_nu7_a02-offset"]
master_directories[47] = ["h04_nu7_a100-offset", "h04_nu7_a50-offset", "h04_nu7_a167-offset", "h04_nu7_a05-offset"]
master_directories[86] = ["h08_nu6_a167-offset", "h08_nu6_a05-offset", "h08_nu6_a02-offset"]
master_directories[66] = ["h06_nu6_a50-offset", "h06_nu6_a167-offset", "h06_nu6_a05-offset"]
master_directories[0] = ["hi_res_high_density-2000", "hi_res_high_density-2000-m2500"]
master_directories[871] = ["h08_nu7_a167-offset", "h08_nu7_a05-offset", "h08_nu7_a02-offset", "h08_nu7_a01-offset", "h08_nu7_a25-low_mass-offset"]
master_directories[671] = ["h06_nu7_a50-offset", "h06_nu7_a167-offset", "h06_nu7_a05-offset", "h06_nu7_a02-offset", "h06_nu7_a125-low_mass-offset"]

master_accretion_rates = {}
master_accretion_rates[87] = [0.17, 0.05, 0.02, 0.01]
master_accretion_rates[67] = [0.50, 0.17, 0.05, 0.02]
master_accretion_rates[47] = [1.00, 0.50, 0.17, 0.05]
master_accretion_rates[86] = [0.17, 0.05, 0.02]
master_accretion_rates[66] = [0.50, 0.17, 0.05]
master_accretion_rates[0] = [0, 0]
master_accretion_rates[871] = [0.17, 0.05, 0.02, 0.01, 0.25]
master_accretion_rates[671] = [0.50, 0.17, 0.05, 0.02, 0.125]

master_start_times = {}
master_start_times[87] = [349, 913, 1751, 2875]
master_start_times[87] = [0, 0, 0, 0]
master_start_times[67] = [108, 217, 451, 788]
master_start_times[47] = [59, 70, 104, 223]
master_start_times[86] = [376, 1064, 0]
master_start_times[66] = [116, 247, 677]
master_start_times[0] = [600, 2500]
master_start_times[871] = [349, 913, 1751, 2875, 0]
master_start_times[671] = [108, 217, 451, 788, 0]

master_end_times = {}
master_end_times[87] = [5000, 4745, 9000, 11700]
master_end_times[67] = [2512, 2502, 6918, 7500]
master_end_times[47] = [2097, 1225, 1898, 2918]
master_end_times[86] = [1816, 2590, 0]
master_end_times[66] = [675, 1336, 1607]
master_end_times[0] = [5000, 3000]
master_end_times[871] = [4000, 4745, 9000, 11700, 0]
master_end_times[671] = [2512, 2502, 6918, 7500, 0]

master_frame_ranges = {}
master_frame_ranges[87] = [[300, 8900, 50], [850, 7400, 50], [1750, 9200, 50], [2850, 13000, 50]]
master_frame_ranges[67] = [[50, 3000, 25], [125, 3500, 25], [300, 7000, 25], [700, 8500, 25]]
master_frame_ranges[47] = [[100, 3000, 25], [100, 3000, 25], [100, 3000, 25], [200, 3000, 25]]
master_frame_ranges[86] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]] 
master_frame_ranges[66] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]] 
master_frame_ranges[0] = [[0, 5000, 25], [2500, 3000, 25]]
master_frame_ranges[871] = [[0, 8000, 200], [0, 7000, 200], [0, 9000, 200], [0, 11700, 200], [0, 2500, 25]]
master_frame_ranges[671] = [[0, 3000, 25], [0, 3000, 25], [0, 7000, 25], [0, 8500, 25], [0, 4000, 25]]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('choice', type = int,
                         help = 'choice of directories')
    parser.add_argument('--dir', dest = "save_directory", default = "rossbyWaveComparisonOverTime",
                         help = 'save directory (default: rossbyWaveComparisonOverTime)')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

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

    parser.add_argument('--negative', dest = "negative", action = 'store_true', default = False,
                         help = 'add negative mass (default: do not)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 22,
                         help = 'fontsize of plot annotations (default: 22)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'fontsize of plot annotations (default: 3)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

directories = master_directories[args.choice]
accretion_rates = master_accretion_rates[args.choice]
start_times = master_start_times[args.choice]
end_times = master_end_times[args.choice]
frame_ranges = master_frame_ranges[args.choice]

### Get Fargo Parameters ###
p = Parameters(directory = "../" + directories[0])

num_rad = p.ny; num_theta = p.nx
r_min = p.ymin; r_max = p.ymax

surface_density_zero = p.sigma0

taper_time = p.masstaper

scale_height = p.aspectratio
viscosity = p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters(directory = "../" + directories[0])

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

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Number of Cores 
num_cores = args.num_cores

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

negative = args.negative

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Data ###

def get_rossby_criteria(args_here):
    # Unwrap Args
    i, frame, directory = args_here

    # Data
    normalized_density = (fromfile("../%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta)) / surface_density_zero

    vrad = (fromfile("../%s/gasvy%d.dat" % (directory, frame)).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
    vtheta = (fromfile("../%s/gasvx%d.dat" % (directory, frame)).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
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

    # Alternate: For inner rossby rad, use minimum in vorticity profile
    inner_rossby_rad, _ = az.get_radial_peak(-1.0 * averaged_vorticity, fargo_par)

    #inner_rossby_rad = rad[inner_max_diff_i]
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


###############################################################################

inner_rossby_rad_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
outer_rossby_rad_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
peak_rad_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))

rossby_rad_difference_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
inner_peak_difference_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
outer_peak_difference_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))

inner_rossby_value_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
outer_rossby_value_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))

###############################################################################

##### PLOTTING #####

colors = ['k', 'cornflowerblue', 'darkorange', 'r', 'green']
labelsize = 19
size = 100

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

def make_plot(show = False):
    # Set up figure
    if args.choice > 0:
        fig = plot.figure(figsize = (7, 6), dpi = dpi)
    else:
        fig = plot.figure(figsize = (7, 2), dpi = dpi)
    ax = fig.add_subplot(111)

    # Iterate
    #max_gap_depth = 0
    for i, directory in enumerate(directories):
        # Frame Range 
        frame_range = util.get_frame_range(frame_ranges[i])

        # Label
        if args.choice > 0:
            scale_height = float(directories[0].split("_")[0][1:]) / 100.0
            log_viscosity = float(directories[0].split("_")[1][2:]) - 2.0
        else:
            scale_height = 0.06
            log_viscosity = 5.0
        accretion_rate = accretion_rates[i]

        start_time = start_times[i]
        end_time = end_times[i]

        #label = r"$h =$ $%.02f$, $\alpha_\mathrm{visc} = 3 \times 10^{-%d}$, A = %.02f" % (scale_height, log_viscosity, accretion_rate)
        if args.choice > 0:
            label = r"$A = %.02f$" % (accretion_rate)
        else:
            labels = ["Default", "Restart"]
            label = labels[i]

        # Data
        #gap_depth_over_time = np.zeros(len(frame_range))

        #for i, frame in enumerate(frame_range):
        #    get_min((i, frame))

        pool_args = [(j, frame, directory) for j, frame in enumerate(frame_range)]

        p = Pool(num_cores)
        p.map(get_rossby_criteria, pool_args)
        p.terminate()

        #if np.max(gap_depth_over_time) > max_gap_depth:
        #    max_gap_depth = np.max(gap_depth_over_time)

        num_frames = len(frame_range)
        this_inner_peak_difference_over_time = np.array(inner_peak_difference_over_time[:num_frames])

        #if i == 4:
        #    this_gap_depth_over_time = this_gap_depth_over_time * (0.3) # low-mass case

        ### Plot ###
        # Basic
        x = frame_range
        y = this_inner_peak_difference_over_time
        result = plot.plot(x, y, c = colors[i], linewidth = linewidth + 1, zorder = 99, label = label, alpha = 0.65)

        # Vortex Lifetime
        if start_time > 1e6:
            # Set to > 0 to use this
            start_time_i = az.my_searchsorted(x, start_time)
            end_time_i = az.my_searchsorted(x, end_time)

            result = plot.plot(x[start_time_i:end_time_i], y[start_time_i:end_time_i], c = colors[i], linewidth = linewidth + 3, zorder = 99)

            plot.scatter(x[start_time_i], y[start_time_i], c = colors[i], s = 150, marker = "o", zorder = 120)
            if args.choice > 0:
                plot.scatter(x[end_time_i], y[end_time_i], c = colors[i], s = 175, marker = "H", zorder = 120)

    # Reference Line
    plot.plot([0, frame_ranges[3][1]], [-0.15, -0.15], c = 'k', linewidth = 1) # Ro = -0.15 is critical value

    if args.choice > 0:
        plot.legend(loc = "lower right", fontsize = fontsize - 4)
    else:
        plot.legend(loc = "upper left", fontsize = fontsize - 4)

    # Axes
    if args.choice > 0:
        x_max = frame_ranges[3][1]
        plot.xlim(0, x_max)
    else:
        x_max = frame_ranges[0][1]
        plot.xlim(0, x_max)

    plot.ylim(0, 0.3)

    unit = "planet orbits"
    plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel("Radii", fontsize = fontsize)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    title = r"$h = %.02f$          $\alpha \approx %s \times 10^{%d}$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2)
    #title = r"$h = %.02f$          $\alpha_\mathrm{disk} = 3 \times 10^{-%d}$" % (scale_height, log_viscosity)
    if args.choice > 0:
        plot.title("%s" % (title), y = 1.015, fontsize = fontsize + 2)

    #top_text = "Incompressible"; top_y = -0.02
    #bottom_text = "Compressible"; bottom_y = -0.38
    #plot.text(0.9 * x_max, top_y, top_text, horizontalalignment = 'right', verticalalignment = 'top', fontsize = fontsize - 3)
    #plot.text(0.1 * x_max, bottom_y, bottom_text, horizontalalignment = 'left', verticalalignment = 'bottom', fontsize = fontsize - 3)

    #title = readTitle()

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/rossbyWaves_choice%d.png" % (save_directory, args.choice)
    else:
        save_fn = "%s/v%04d_rossbyWaves_choice%d.png" % (save_directory, version, args.choice)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
