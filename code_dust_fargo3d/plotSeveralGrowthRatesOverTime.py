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
from scipy.ndimage import filters as ff

import matplotlib
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import NullFormatter

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
master_start_times[67] = [20, 50, 200, 575]
master_start_times[47] = [59, 70, 104, 223]
master_start_times[86] = [376, 1064, 0]
master_start_times[66] = [116, 247, 677]
master_start_times[0] = [600, 2500]
master_start_times[871] = [150, 720, 1600, 2700, 550]
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
master_frame_ranges[87] = [[300, 8900, 50], [850, 8500, 50], [1750, 9200, 50], [2850, 13000, 50]]
master_frame_ranges[67] = [[0, 3000, 25], [0, 3000, 25], [0, 7000, 25], [0, 8500, 25]]
master_frame_ranges[47] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25], [0, 3000, 25]]
master_frame_ranges[86] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]] 
master_frame_ranges[66] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]] 
master_frame_ranges[0] = [[0, 5000, 25], [2500, 3000, 25]]
master_frame_ranges[871] = [[0, 3900, 25], [0, 4700, 25], [0, 8900, 25], [0, 12500, 25], [0, 7500, 25]]
master_frame_ranges[671] = [[0, 3000, 25], [0, 3000, 25], [0, 7000, 25], [0, 8500, 25], [0, 4000, 25]]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('choice', type = int,
                         help = 'choice of directories')
    parser.add_argument('--dir', dest = "save_directory", default = "growthRates",
                         help = 'save directory (default: vorticity)')
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

### Helper Functions ###
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

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

### Data ###

def get_contrasts(args_here):
    # Unwrap Args
    i, frame, directory = args_here

    if frame == 8800:
        frame = 8801 # Remove problem frame

    # Get Data
    density = fromfile("../%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta) / surface_density_zero
    density, shift_c = shift_density(density, fargo_par, reference_density = density)
    azimuthal_profile = az.get_mean_azimuthal_profile(density, fargo_par, sliver_width = 1.5, end = 1.85)

    if "low_mass" in directory:
        azimuthal_profile /= (0.3) # low-mass case

    maxima_over_time[i] = np.percentile(azimuthal_profile, 90)
    minima_over_time[i] = np.percentile(azimuthal_profile, 5)
    contrasts_over_time[i] = maxima_over_time[i] / minima_over_time[i]
    differences_over_time[i] = maxima_over_time[i] - minima_over_time[i]

    print i, frame, differences_over_time[i], maxima_over_time[i], minima_over_time[i], contrasts_over_time[i]

###############################################################################

maxima_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
minima_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
contrasts_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
differences_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))

###############################################################################

##### PLOTTING #####

colors = ['k', 'cornflowerblue', 'darkorange', 'r', 'purple']
labelsize = 19
size = 100
alpha = 0.8

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

def make_plot(show = False):
    # Set up figure
    fig, (ax1, ax2, ax3) = plot.subplots(3, 1, figsize = (6, 12), gridspec_kw={'height_ratios': [4, 6, 5]})

    # Iterate
    for i, directory in enumerate(directories):
        # Frame Range 
        frame_range = util.get_frame_range(frame_ranges[i])
        dt = (frame_range[1] - frame_range[0]) * (2.0 * np.pi) # for growth rate calculation

        start_time = start_times[i]
        start_time_i = az.my_searchsorted(frame_range, start_time)

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
            if i == 4:
                master_label = r"$0.3$ $\Sigma_0$"
            else:
                master_label = r"$A = %.02f$" % (accretion_rate)
        else:
            labels = ["Default", "Restart"]
            label = labels[i]

        # Data
        #gap_depth_over_time = np.zeros(len(frame_range))

        #for i, frame in enumerate(frame_range):
        #    get_min((i, frame))

        pool_args = [(j, frame, directory) for j, frame in enumerate(frame_range)]

        p = Pool(num_cores)
        p.map(get_contrasts, pool_args)
        p.terminate()

        num_frames = len(frame_range)
        this_maxima_over_time = np.array(maxima_over_time[:num_frames])
        this_minima_over_time = np.array(minima_over_time[:num_frames])
        this_contrasts_over_time = np.array(contrasts_over_time[:num_frames])
        this_differences_over_time = np.array(differences_over_time[:num_frames])

        this_smoothed_differences_over_time = smooth(this_differences_over_time, 5)
        this_growth_rates_over_time = np.diff(np.log(this_smoothed_differences_over_time)) / dt

        ##### Top Plot #####

        # Plot
        x = frame_range
        y1 = this_contrasts_over_time
        p1, = ax1.plot(x[start_time_i:], y1[start_time_i:], c = colors[i], linewidth = linewidth, zorder = 99 - i)

        # Axes
        ax1.set_yscale('log')

        ax1.yaxis.set_major_formatter(ScalarFormatter())
        ax1.yaxis.set_minor_formatter(NullFormatter())

        if i == 3:
            ax1.set_xlim(x[0], x[-1])

        if scale_height == 0.08:
            ax1.set_ylim(1, 15)
            ax1.set_yticks([1, 3, 10])
        else:
            ax1.set_ylim(1, 7)
            ax1.set_yticks([1, 3, 7])

        # Annotate
        #ax1.set_xlabel("", fontsize = fontsize)
        ax1.set_ylabel("Contrast", fontsize = fontsize)

        alpha_coefficent = "3"
        if scale_height == 0.08:
            alpha_coefficent = "1.5"
        elif scale_height == 0.04:
            alpha_coefficent = "6"

        title1 = r"$h = %.2f$    $\alpha = %s \times 10^{%d}$" % (scale_height, alpha_coefficent, int(round(np.log(viscosity) / np.log(10), 0)) + 2)
        #title1 = r"$A = %.2f$" % (accretion)
        ax1.set_title("%s" % (title1), y = 1.035, fontsize = fontsize + 2)

        ##### Middle Plot #####

        label1 = ""; label2 = ""
        if i == 0:
            label1 = r"$\Sigma_\mathrm{max}$"
            label2 = r"$\Sigma_\mathrm{min}$"

        y2 = this_maxima_over_time
        p2, = ax2.plot(x[start_time_i:], y2[start_time_i:], c = colors[i], linewidth = linewidth, label = label1, zorder = 99 - i)
        y3 = this_minima_over_time
        p3, = ax2.plot(x[start_time_i:], y3[start_time_i:], c = colors[i], linewidth = linewidth - 2, label = label2, zorder = 90 - i)

        # Axes
        if i == 3:
            ax2.set_xlim(x[0], x[-1])

        if scale_height == 0.08:
            ax2.set_ylim(0, 2.8)
            ax2.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
        else:
            ax2.set_ylim(0, 1.8)
            ax2.set_yticks([0, 0.5, 1, 1.5])

        # Annotate
        ax2.set_ylabel(r"$\Sigma$ $/$ $\Sigma_0$", fontsize = fontsize)
        ax2.legend(loc = "upper right", fontsize = fontsize - 5)

        ##### Bottom Plot #####

        # Plot
        y4 = this_growth_rates_over_time
        p4, = ax3.plot(x[start_time_i:-5], y4[start_time_i:-4], c = colors[i], linewidth = linewidth, alpha = alpha, label = master_label, zorder = 90 - i)

        # Axes
        if i == 3:
            ax3.set_xlim(x[0], x[-1])
        ax3.set_ylim(0.4 * 10**(-5), 0.4 * 10**(-2))
        ax3.set_yscale('log')

        ax3.set_yticks([10**(-5), 10**(-4), 10**(-3)])

        # Annotate
        ax3.set_xlabel("Time (planet orbits)", fontsize = fontsize)
        ax3.set_ylabel(r"Growth Rate", fontsize = fontsize)
        ax3.legend(loc = "upper right", fontsize = fontsize - 5)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/growthRatesOverTime_choice%d.png" % (save_directory, args.choice)
    else:
        save_fn = "%s/v%04d_growthRatesOverTime_choice%d.png" % (save_directory, version, arg.choice)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
