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
import matplotlib.style
#matplotlib.style.use('classic')
from matplotlib.ticker import ScalarFormatter
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
#master_start_times[67] = [108, 217, 451, 788]
master_start_times[67] = [75, 175, 380, 770]
#master_start_times[47] = [59, 70, 104, 223]
master_start_times[47] = [36, 40, 58, 75]
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
master_frame_ranges[87] = [[0, 8000, 50], [0, 7000, 50], [0, 9000, 50], [0, 11700, 50]]
master_frame_ranges[67] = [[50, 3000, 25], [75, 3500, 25], [100, 7000, 25], [200, 8500, 25]]
master_frame_ranges[47] = [[25, 3000, 25], [50, 3000, 25], [50, 3000, 25], [50, 3000, 25]]
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
    parser.add_argument('--dir', dest = "save_directory", default = "pressureBumps",
                         help = 'save directory (default: pressureBumps)')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

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

# Initial Vortices

initial_accretion_numbers_h4 = np.array([1, 2, 3, 4])
initial_times_h4 = np.array([36, 40, 58, 75])
initial_critical_maxima_h4 = np.array([1.120, 1.127, 1.130, 1.117])
initial_density_maxima_h4 = np.array([1.160, 1.160, 1.156, 1.146])

initial_accretion_numbers_h6 = np.array([1, 2, 3, 4])
initial_times_h6 = np.array([75, 175, 380, 770])
initial_critical_maxima_h6 = np.array([1.168, 1.182, 1.200, 1.217])
initial_density_maxima_h6 = np.array([1.220, 1.228, 1.242, 1.257])

# Later Generation Vortices

second_accretion_numbers_h4 = np.array([2, 3, 3, 4, 4, 4, 4])
second_times_h4 = np.array([790, 900, 1220, 1010, 1230, 1530, 1770])
second_critical_maxima_h4 = np.array([1.243, 1.192, 1.204, 1.172, 1.164, 1.174, 1.178])
second_density_maxima_h4 = np.array([1.343, 1.300, 1.310, 1.250, 1.257, 1.264, 1.275])

second_accretion_numbers_h6 = np.array([3, 3, 4, 4, 4])
second_times_h6 = np.array([2370, 3400, 3080, 3550, 3910])
second_critical_maxima_h6 = np.array([1.228, 1.228, 1.250, 1.248, 1.304])
second_density_maxima_h6 = np.array([1.396, 1.414, 1.364, 1.375, 1.382])

# Vortices that reform interior

interior_accretion_numbers_h4 = np.array([4])
interior_times_h4 = np.array([2140])
interior_critical_maxima_h4 = np.array([1.193])
interior_density_maxima_h4 = np.array([1.304])

interior_accretion_numbers_h6 = np.array([3, 3])
interior_times_h6 = np.array([4300, 6000])
interior_critical_maxima_h6 = np.array([1.224, 1.235])
interior_density_maxima_h6 = np.array([1.432, 1.457])

# Vortices that reform too interior

too_interior_accretion_numbers_h4 = np.array([1, 2, 3, 4])
too_interior_times_h4 = np.array([900, 1050, 1700, 2950])
too_interior_critical_maxima_h4 = np.array([1.202, 1.200, 1.189, 1.164])
too_interior_density_maxima_h4 = np.array([1.382, 1.357, 1.325, 1.304])

too_interior_accretion_numbers_h6 = np.array([1, 2])
too_interior_times_h6 = np.array([2480, 1760])
too_interior_critical_maxima_h6 = np.array([1.264, 1.256])
too_interior_density_maxima_h6 = np.array([1.540, 1.465])

###############################################################################

### Helper Functions ###

def range_brace(x_min, x_max, mid=0.75, 
                beta1=50.0, beta2=100.0, height=1, 
                initial_divisions=11, resolution_factor=1.5):
    # Source: http://stackoverflow.com/questions/1289681/drawing-braces-with-pyx
    # determine x0 adaptively values using second derivitive
    # could be replaced with less snazzy:
    #   x0 = np.arange(0, 0.5, .001)
    x0 = np.array(())
    tmpx = np.linspace(0, 0.5, initial_divisions)
    tmp = beta1**2 * (np.exp(beta1*tmpx)) * (1-np.exp(beta1*tmpx)) / np.power((1+np.exp(beta1*tmpx)),3)
    tmp += beta2**2 * (np.exp(beta2*(tmpx-0.5))) * (1-np.exp(beta2*(tmpx-0.5))) / np.power((1+np.exp(beta2*(tmpx-0.5))),3)
    for i in range(0, len(tmpx)-1):
        t = int(np.ceil(resolution_factor*max(np.abs(tmp[i:i+2]))/float(initial_divisions)))
        x0 = np.append(x0, np.linspace(tmpx[i],tmpx[i+1],t))
    x0 = np.sort(np.unique(x0)) # sort and remove dups
    # half brace using sum of two logistic functions
    y0 = mid*2*((1/(1.+np.exp(-1*beta1*x0)))-0.5)
    y0 += (1-mid)*2*(1/(1.+np.exp(-1*beta2*(x0-0.5))))
    # concat and scale x
    x = np.concatenate((x0, 1-x0[::-1])) * float((x_max-x_min)) + x_min
    y = np.concatenate((y0, y0[::-1])) * float(height)
    return (x,y)

def find_min(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 0.9) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 1.25) # look for max radial density before r = 2.3
    min_rad_outer_index = np.argmin(averagedDensity[outer_disk_start : outer_disk_end])

    min_index = outer_disk_start + min_rad_outer_index
    min_rad = rad[min_index]
    min_density = averagedDensity[min_index]

    return min_density

### Data ###

def get_min(args_here):
    # Unwrap Args
    i, frame, directory = args_here

    # Get Data
    density = fromfile("../%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta)
    avg_density = np.average(density, axis = 1)
    normalized_density = avg_density / surface_density_zero
    
    # Get Minima
    radial_peak_a, _ = az.get_radial_peak(avg_density, fargo_par, end = 1.6)

    # Print Update
    print "%d: %.3f" % (frame, radial_peak_a)

    # Store Data
    radial_peak_over_time[i] = radial_peak_a

###############################################################################

radial_peak_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))

###############################################################################

##### PLOTTING #####

colors = np.array(['k', 'cornflowerblue', 'darkorange', 'r', 'green'])
labelsize = 19
alpha = 1.0
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
    max_gap_depth = 0
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
        p.map(get_min, pool_args)
        p.terminate()

        #if np.max(gap_depth_over_time) > max_gap_depth:
        #    max_gap_depth = np.max(radial_peak_over_time)

        num_frames = len(frame_range)
        this_radial_peak_over_time = np.array(radial_peak_over_time[:num_frames])

        ### Plot ###
        # Basic
        x = frame_range
        y = this_radial_peak_over_time
        result = plot.plot(x, y, c = colors[i], linewidth = linewidth - 1, zorder = 99, label = label)

        # Vortex Lifetime
        if start_time > 0:
            start_time_i = az.my_searchsorted(x, start_time)
            end_time_i = az.my_searchsorted(x, end_time)

            result = plot.plot(x[start_time_i:end_time_i], y[start_time_i:end_time_i], c = colors[i], linewidth = linewidth + 3, zorder = 99)

            #plot.scatter(x[start_time_i], y[start_time_i], c = colors[i], s = 150, marker = "o", zorder = 120)
            if args.choice > 0:
                plot.scatter(x[end_time_i], y[end_time_i], c = colors[i], s = 175, marker = "H", zorder = 120)

    # Scatter critical points
    if scale_height == 0.04:
        i_h4 = plot.scatter(initial_times_h4, initial_critical_maxima_h4, s = 120, c = colors[initial_accretion_numbers_h4 - 1], zorder = 100, alpha = alpha)
        s_h4 = plot.scatter(second_times_h4, second_critical_maxima_h4, s = 170, c = colors[second_accretion_numbers_h4 - 1], zorder = 100, alpha = alpha, marker = "*")
        d_h4 = plot.scatter(interior_times_h4, interior_critical_maxima_h4, s = 120, c = colors[interior_accretion_numbers_h4 - 1], zorder = 100, alpha = alpha, marker = "D")
        t_h4 = plot.scatter(too_interior_times_h4, too_interior_critical_maxima_h4, s = 120, c = colors[too_interior_accretion_numbers_h4 - 1], zorder = 100, alpha = alpha, marker = "x")
    elif scale_height == 0.06:
        i_h6 = plot.scatter(initial_times_h6, initial_critical_maxima_h6, s = 120, c = colors[initial_accretion_numbers_h6 - 1], zorder = 100, alpha = alpha)
        s_h4 = plot.scatter(second_times_h6, second_critical_maxima_h6, s = 170, c = colors[second_accretion_numbers_h6 - 1], zorder = 100, alpha = alpha, marker = "*")
        d_h4 = plot.scatter(interior_times_h6, interior_critical_maxima_h6, s = 120, c = colors[interior_accretion_numbers_h6 - 1], zorder = 100, alpha = alpha, marker = "D")
        t_h4 = plot.scatter(too_interior_times_h6, too_interior_critical_maxima_h6, s = 120, c = colors[too_interior_accretion_numbers_h6 - 1], zorder = 100, alpha = alpha, marker = "x")

    if scale_height == 0.04:
        plot.legend(loc = "upper left", fontsize = fontsize - 4)

        # Pressure Bump label
        x1 = 1200; x2 = 2800; x_center = x1 + (x2 - x1) / 2.0
        y_base = 1.45; dy = 0.08

        plot.text(x_center, y_base + dy, r"$r_\mathrm{pressure}$", horizontalalignment = 'center', fontsize = fontsize)
        top_brace_x, top_brace_y = range_brace(x1, x2, height = 0.06)
        plot.plot(top_brace_x, y_base + top_brace_y, c = "k", linewidth = 2)

        # Critical Bump label
        x1 = 800; x2 = 2200; x_center = x1 + (x2 - x1) / 2.0
        y_base = 1.14; dy = 0.10

        plot.text(x_center, y_base - dy, r"$r_\mathrm{crit}$", horizontalalignment = 'center', fontsize = fontsize)
        bottom_brace_x, bottom_brace_y = range_brace(x1, x2, height = 0.06)
        plot.plot(bottom_brace_x, y_base - bottom_brace_y, c = "k", linewidth = 2)

    elif scale_height == 0.06:
        plot.legend(loc = "lower right", fontsize = fontsize - 4)

        # Pressure Bump label
        x1 = 4000; x2 = 6000; x_center = x1 + (x2 - x1) / 2.0
        y_base = 1.48; dy = 0.08

        plot.text(x_center, y_base + dy, r"$r_\mathrm{pressure}$", horizontalalignment = 'center', fontsize = fontsize)
        top_brace_x, top_brace_y = range_brace(x1, x2, height = 0.06)
        plot.plot(top_brace_x, y_base + top_brace_y, c = "k", linewidth = 2)

        # Critical Bump label
        x1 = 1800; x2 = 4200; x_center = x1 + (x2 - x1) / 2.0
        y_base = 1.18; dy = 0.10

        plot.text(x_center, y_base - dy, r"$r_\mathrm{crit}$", horizontalalignment = 'center', fontsize = fontsize)
        bottom_brace_x, bottom_brace_y = range_brace(x1, x2, height = 0.06)
        plot.plot(bottom_brace_x, y_base - bottom_brace_y, c = "k", linewidth = 2)


    # Axes
    if args.choice > 0:
        plot.xlim(0, frame_range[-1])
    else:
        plot.xlim(0, frame_ranges[0][1])

    start_y = 1.0; end_y = 1.6
    plot.ylim(start_y, end_y)

    #title = readTitle()

    unit_x = "planet orbits"; unit_y = "r_\mathrm{p}"
    plot.xlabel(r"Time [%s]" % unit_x, fontsize = fontsize)
    plot.ylabel(r"$r$ [$%s$]" % unit_y, fontsize = fontsize)

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
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    #title = readTitle()

    # Set twin axis
    twin = ax.twinx()
    twin.set_ylim(0, (end_y - start_y) / scale_height)
    twin.set_ylabel(r"$(r - r_\mathrm{p})$ $/$ $h$", fontsize = fontsize, rotation = 270, labelpad = 30)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/pressureBumps_choice%d.png" % (save_directory, args.choice)
    else:
        save_fn = "%s/v%04d_pressureBumps_choice%d.png" % (save_directory, version, arg.choice)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
