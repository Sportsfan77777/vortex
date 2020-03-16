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

master_directories = {}
master_directories[87] = ["h08_nu7_a167-offset", "h08_nu7_a05-offset", "h08_nu7_a02-offset", "h08_nu7_a01-offset"]
master_directories[67] = ["h06_nu7_a50-offset", "h06_nu7_a167-offset", "h06_nu7_a05-offset", "h06_nu7_a02-offset"]
master_directories[47] = ["h04_nu7_a100-offset", "h04_nu7_a50-offset", "h04_nu7_a167-offset", "h04_nu7_a05-offset"]
master_directories[86] = ["h08_nu6_a167-offset", "h08_nu6_a05-offset", "h08_nu6_a02-offset"]
master_directories[66] = ["h06_nu6_a50-offset", "h06_nu6_a167-offset", "h06_nu6_a05-offset"]
master_directories[0] = ["hi_res_high_density-2000"]
master_directories[871] = ["h08_nu7_a167-offset", "h08_nu7_a05-offset", "h08_nu7_a02-offset", "h08_nu7_a01-offset", "h08_nu7_a25-low_mass-offset"]
master_directories[671] = ["h06_nu7_a50-offset", "h06_nu7_a167-offset", "h06_nu7_a05-offset", "h06_nu7_a02-offset", "h06_nu7_a125-low_mass-offset"]

master_accretion_rates = {}
master_accretion_rates[87] = [0.17, 0.05, 0.02, 0.01]
master_accretion_rates[67] = [0.50, 0.17, 0.05, 0.02]
master_accretion_rates[47] = [1.00, 0.50, 0.17, 0.05]
master_accretion_rates[86] = [0.17, 0.05, 0.02]
master_accretion_rates[66] = [0.50, 0.17, 0.05]
master_accretion_rates[0] = [0]
master_accretion_rates[871] = [0.17, 0.05, 0.02, 0.01, 0.25]
master_accretion_rates[671] = [0.50, 0.17, 0.05, 0.02, 0.125]

master_start_times = {}
master_start_times[87] = [349, 913, 1751, 2875]
master_start_times[67] = [108, 217, 451, 788]
master_start_times[47] = [59, 70, 104, 223]
master_start_times[86] = [376, 1064, 0]
master_start_times[66] = [116, 247, 677]
master_start_times[0] = [600]
master_start_times[871] = [349, 913, 1751, 2875, 0]
master_start_times[671] = [108, 217, 451, 788, 0]

master_end_times = {}
master_end_times[87] = [4000, 4745, 6790, 10700]
master_end_times[67] = [2512, 2502, 6918, 7500]
master_end_times[47] = [2097, 1225, 1898, 2918]
master_end_times[86] = [1816, 2590, 0]
master_end_times[66] = [675, 1336, 1607]
master_end_times[0] = [4000]
master_end_times[871] = [4000, 4745, 6790, 10700, 0]
master_end_times[671] = [2512, 2502, 6918, 7500, 0]

master_frame_ranges = {}
master_frame_ranges[87] = [[400, 3500, 250], [1000, 4200, 250], [1800, 8300, 250], [3000, 11700, 250]]
master_frame_ranges[67] = [[200, 3000, 250], [300, 3000, 250], [500, 6800, 250], [800, 7500, 250]]
master_frame_ranges[47] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25], [0, 3000, 25]]
master_frame_ranges[86] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]] 
master_frame_ranges[66] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]]
master_frame_ranges[0] = [[0, 4000, 25]]
master_frame_ranges[871] = [[0, 8000, 200], [0, 7000, 200], [0, 7000, 200], [0, 11700, 200], [0, 2500, 25]]
master_frame_ranges[671] = [[0, 3000, 25], [0, 3000, 25], [0, 7000, 25], [0, 8500, 25], [0, 4000, 25]]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('choice', type = int,
                         help = 'choice of directories')
    parser.add_argument('--dir', dest = "save_directory", default = "stress",
                         help = 'save directory (default: stress)')
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

### Helper Functions ###

def get_reynolds_stress(args_here):
    # Unwrap Args
    i, frame, directory = args_here

    # Data
    density = fromfile("../%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta) / surface_density_zero
    averagedDensity = np.average(density, axis = 1)

    radial_velocity = fromfile("../%s/gasvy%d.dat" % (directory, frame)).reshape(num_rad, num_theta)
    azimuthal_velocity = fromfile("../%s/gasvx%d.dat" % (directory, frame)).reshape(num_rad, num_theta)

    #keplerian_velocity = rad * (np.power(rad, -1.5) - 1) # in rotating frame, v_k = r * (r^-1.5 - r_p^-1.5)
    #sub_keplerian_velocity = keplerian_velocity - 0.5 * np.power(scale_height, 2)
    #azimuthal_velocity -= sub_keplerian_velocity[:, None]

    radial_velocity -= np.average(radial_velocity, axis = 1)[:, None]
    azimuthal_velocity -= np.average(azimuthal_velocity, axis = 1)[:, None]

    sound_speed = scale_height * np.power(rad, -1.5)

    stress = np.multiply(radial_velocity, azimuthal_velocity)
    averagedStress = np.abs(np.average(stress, axis = 1) / np.power(sound_speed, 2))

    # Zoom in on region around vortex
    vortex_left_rad, vortex_right_rad = az.get_radial_bounds(density, fargo_par, threshold = 0.9)
    vortex_left_i = np.searchsorted(rad, vortex_left_rad)
    vortex_right_i = np.searchsorted(rad, vortex_right_rad)
    zoom_averagedStress = averagedStress[vortex_left_i : vortex_right_i]

    median_averagedStress = np.median(zoom_averagedStress)
    max_averagedStress = np.max(zoom_averagedStress)

    # Print Update
    print "%d: %.3f, %.3f" % (frame, median_averagedStress, max_averagedStress)

    # Store Data
    stress_over_time[i] = median_averagedStress

###############################################################################

stress_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))

###############################################################################

##### PLOTTING #####

colors = ['k', 'cornflowerblue', 'darkorange', 'r', 'green']
labelsize = 19
size = 100

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Iterate
    max_stress = 0
    for i, directory in enumerate(directories):
        # Frame Range 
        frame_range = util.get_frame_range(frame_ranges[i])

        # Label
        scale_height = float(directories[0].split("_")[0][1:]) / 100.0
        log_viscosity = float(directories[0].split("_")[1][2:]) - 2.0
        accretion_rate = accretion_rates[i]

        start_time = start_times[i]
        end_time = end_times[i]

        #label = r"$h =$ $%.02f$, $\alpha_\mathrm{visc} = 3 \times 10^{-%d}$, A = %.02f" % (scale_height, log_viscosity, accretion_rate)
        label = r"$A = %.02f$" % (accretion_rate)

        # Data
        #get_reynolds_stress = np.zeros(len(frame_range))

        #for i, frame in enumerate(frame_range):
        #    get_reynolds_stress((i, frame, directory))

        pool_args = [(j, frame, directory) for j, frame in enumerate(frame_range)]

        p = Pool(num_cores)
        p.map(get_reynolds_stress, pool_args)
        p.terminate()

        if np.max(stress_over_time) > max_stress:
            max_stress = np.max(stress_over_time)

        num_frames = len(frame_range)
        this_stress_over_time = np.array(stress_over_time[:num_frames])

        if i == 4:
            this_stress_over_time = this_stress_over_time * (0.3) # low-mass case

        ### Plot ###
        # Basic
        x = frame_range
        y = this_stress_over_time
        result = plot.plot(x, y, c = colors[i], linewidth = linewidth - 1, zorder = 99, label = label)

        # Vortex Lifetime
        if start_time > 0:
            start_time_i = az.my_searchsorted(x, start_time)
            end_time_i = az.my_searchsorted(x, end_time)

            result = plot.plot(x[start_time_i:end_time_i], y[start_time_i:end_time_i], c = colors[i], linewidth = linewidth + 3, zorder = 99)

            plot.scatter(x[start_time_i], y[start_time_i], c = colors[i], s = 150, marker = "o", zorder = 120)
            plot.scatter(x[end_time_i], y[end_time_i], c = colors[i], s = 175, marker = "H", zorder = 120)

    plot.legend(loc = "upper right", fontsize = fontsize - 4)

    # Axes
    plot.xlim(0, frame_range[-1])
    plot.ylim(10**(-5), 3 * 10**(-1))

    plot.yscale('log')

    #title = readTitle()

    unit = "planet orbits"
    plot.xlabel(r"Time [%s]" % unit, fontsize = fontsize)
    plot.ylabel(r"Reynolds Stress ($<\Delta v_r \Delta v_{\phi}> / c_s^2$)", fontsize = fontsize)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    title = r"$h = %.02f$          $\alpha = %s \times 10^{%d}$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2)
    #title = r"$h = %.02f$          $\alpha_\mathrm{disk} = 3 \times 10^{-%d}$" % (scale_height, log_viscosity)
    plot.title("%s" % (title), y = 1.015, fontsize = fontsize + 2)
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    #title = readTitle()

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/stress_choice%d.png" % (save_directory, args.choice)
    else:
        save_fn = "%s/v%04d_stress_choice%d.png" % (save_directory, version, arg.choice)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
