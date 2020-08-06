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

problem_frames = np.array([3000, 3125, 3200, 3250, 3550, 7425]) # h = 0.08, A = 0.167

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
master_frame_ranges[87] = [[400, 10000, 25], [950, 7400, 25], [1750, 9200, 25], [2900, 12350, 25]]
master_frame_ranges[67] = [[0, 3000, 25], [0, 3000, 25], [0, 7000, 25], [0, 8500, 25]]
master_frame_ranges[47] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25], [0, 3000, 25]]
master_frame_ranges[86] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]] 
master_frame_ranges[66] = [[0, 3000, 25], [0, 3000, 25], [0, 3000, 25]] 
master_frame_ranges[0] = [[0, 5000, 25], [2500, 3000, 25]]
master_frame_ranges[871] = [[400, 10000, 25], [950, 7400, 25], [1750, 9200, 25], [2900, 12350, 25], [600, 7500, 25]]
master_frame_ranges[671] = [[0, 3000, 25], [0, 3000, 25], [0, 7000, 25], [0, 8500, 25], [0, 4000, 25]]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('choice', type = int,
                         help = 'choice of directories')
    parser.add_argument('--dir', dest = "save_directory", default = "extents",
                         help = 'save directory (default: extents)')
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

    parser.add_argument('-r', dest = "check_rossby", type = int, default = 1000000,
                         help = 'frame at which you start using the Rossby number for measuring everything (default: infinity)')
    parser.add_argument('-e', dest = "extreme_cutoff", type = int, default = 1000000,
                         help = 'frame at which you start using the extreme Rossby number cutoff (default: infinity)')
    parser.add_argument('-i', dest = "include_aspect", action = 'store_true', default = False,
                         help = 'include aspect ratio (default: do not)')
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

check_rossby = args.check_rossby
extreme_cutoff = args.extreme_cutoff
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

### Data ###

check_rossby_frames = [2000, 2700, 8000, 12100, 100000]

def get_extents(args_here):
    # Unwrap Args
    i, frame, directory, check_rossby = args_here

    if frame in problem_frames:
        frame += 3 # switch to an adjacent frame

    # Get Data
    density = fromfile("../%s/gasdens%d.dat" % (directory, frame)).reshape(num_rad, num_theta) / surface_density_zero
    avg_density = np.average(density, axis = 1)
    peak_rad, peak_density = az.get_radial_peak(avg_density, fargo_par)

    normal = True
    if frame > check_rossby:
        vrad = (fromfile("../%s/gasvy%d.dat" % (directory, frame)).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
        vtheta = (fromfile("../%s/gasvx%d.dat" % (directory, frame)).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
        vorticity = utilVorticity.velocity_curl(vrad, vtheta, rad, theta, rossby = True, residual = True)

        # Find minimum
        if accretion > 0.015:
            start_rad = min([peak_rad - 0.05, 1.5])
            start_rad_i = np.searchsorted(rad, start_rad) # Is this necessary?
            end_rad_i = np.searchsorted(rad, 2.5)
        else:
            start_rad_i = np.searchsorted(rad, 1.0) # Is this necessary?
            end_rad_i = np.searchsorted(rad, 1.8)
        zoom_vorticity = vorticity[start_rad_i : end_rad_i]

        min_rossby_number = np.percentile(zoom_vorticity, 0.25)
        if min_rossby_number < -0.15:
            normal = False # Compressible regime from Surville+ 15

    if normal:
        azimuthal_extent = az.get_extent(density, fargo_par, threshold = 0.85) # Use 0.9 for h = 0.08 (Add as a parameter)
        radial_extent, radial_peak = az.get_radial_extent(density, fargo_par, threshold = 0.85)
        radial_peak_a, _ = az.get_radial_peak(avg_density, fargo_par)

        azimuthal_extent_over_time[i] = azimuthal_extent * (180.0 / np.pi)
        radial_extent_over_time[i] = radial_extent / scale_height
        radial_peak_over_time[i] = radial_peak_a # radial_peak
        #radial_peak_over_time_a[i] = radial_peak_a
    else:
        # Shift everything
        density, vorticity, shift_c = shift_density(density, vorticity, fargo_par, reference_density = density)

        # Locate minimum
        if accretion > 0.015:
            start_rad = min([peak_rad - 0.05, 1.5])
            start_rad_i = np.searchsorted(rad, start_rad) # Is this necessary?
            end_rad_i = np.searchsorted(rad, 2.5)
        else:
            start_rad_i = np.searchsorted(rad, 1.0) # Is this necessary?
            end_rad_i = np.searchsorted(rad, 1.8)
        zoom_vorticity = vorticity[start_rad_i : end_rad_i]

        min_rossby_number = np.percentile(zoom_vorticity, 0.25)
        abs_zoom_vorticity = np.abs(zoom_vorticity - min_rossby_number)
        minimum_location = np.argmin(abs_zoom_vorticity)

        rad_min_i, theta_min_i = np.unravel_index(minimum_location, np.shape(zoom_vorticity))

        # Locate radial and azimuthal center
        left_side = zoom_vorticity[rad_min_i, :theta_min_i]
        right_side = zoom_vorticity[rad_min_i, theta_min_i:]
        front_side = zoom_vorticity[:rad_min_i, theta_min_i]
        back_side = zoom_vorticity[rad_min_i:, theta_min_i]

        if frame < extreme_cutoff:
            cutoff = -0.04
        else:
            cutoff = -0.12 # Extreme! (neglects "vortex" that develops around the minimum)

        left_i = theta_min_i - az.my_searchsorted(left_side[::-1], cutoff) # at location of minimum
        right_i = theta_min_i + az.my_searchsorted(right_side, cutoff)
        front_i = rad_min_i - az.my_searchsorted(front_side[::-1], cutoff)
        back_i = rad_min_i + az.my_searchsorted(back_side, cutoff)

        radial_center_i = start_rad_i + int((front_i + back_i) / 2.0)
        azimuthal_center_i = int((left_i + right_i) / 2.0)

        radial_center = (rad[start_rad_i + front_i] + rad[start_rad_i + back_i]) / 2.0
        azimuthal_center = ((theta[left_i] + theta[right_i]) / 2.0) * (180.0 / np.pi)

        print i, frame, rad[start_rad_i + rad_min_i], theta[theta_min_i] * (180.0 / np.pi), "Minimum Rossby Number"
        print i, frame, rad[start_rad_i + front_i], radial_center, rad[start_rad_i + back_i], "Radial: Left, Center, Right"
        print i, frame, theta[left_i] * (180.0 / np.pi), azimuthal_center, theta[right_i] * (180.0 / np.pi), "Azimuthal: Left, Center, Right"

        # Measure radial and azimuthal extents
        left_side = vorticity[radial_center_i, :azimuthal_center_i]
        right_side = vorticity[radial_center_i, azimuthal_center_i:]
        front_side = vorticity[:radial_center_i, azimuthal_center_i]
        back_side = vorticity[radial_center_i:, azimuthal_center_i]

        left_i = azimuthal_center_i - az.my_searchsorted(left_side[::-1], cutoff) # relative to center
        right_i = azimuthal_center_i + az.my_searchsorted(right_side, cutoff)
        front_i = radial_center_i - az.my_searchsorted(front_side[::-1], cutoff)
        back_i = radial_center_i + az.my_searchsorted(back_side, cutoff)

        radial_peak_over_time[i] = radial_center
        radial_extent_over_time[i] = (rad[back_i] - rad[front_i]) / scale_height
        azimuthal_extent_over_time[i] = theta[right_i - left_i] * (180.0 / np.pi)

        print i, frame, rad[front_i], radial_center, rad[back_i], "Final Radial: Left, Center, Right"
        print i, frame, theta[left_i] * (180.0 / np.pi), azimuthal_center, theta[right_i] * (180.0 / np.pi), "Final Azimuthal: Left, Center, Right"

    #contrasts_over_time[i] = az.get_contrast(density, fargo_par)

    print i, frame, azimuthal_extent_over_time[i], radial_extent_over_time[i], radial_peak_over_time[i]
    #print i, frame, azimuthal_extent_over_time[i], radial_extent_over_time[i], radial_peak_over_time[i], radial_peak_over_time_a[i], contrasts_over_time[i]


###############################################################################

azimuthal_extent_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
radial_extent_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))
radial_peak_over_time = mp_array("d", 10 * len(util.get_frame_range(frame_ranges[0])))

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
    fig, (ax1, ax2, ax3) = plot.subplots(3, 1, figsize = (6, 12), gridspec_kw={'height_ratios': [1, 1, 1]})

    # Iterate
    for i, directory in enumerate(directories):
        # Frame Range 
        frame_range = util.get_frame_range(frame_ranges[i])
        dt = (frame_range[1] - frame_range[0]) * (2.0 * np.pi) # for growth rate calculation

        start_time = start_times[i]
        start_time_i = az.my_searchsorted(frame_range, start_time)

        check_rossby = check_rossby_frames[i]

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
        #for j, frame in enumerate(frame_range):
        #    get_extents((j, frame, directory))

        pool_args = [(j, frame, directory, check_rossby) for j, frame in enumerate(frame_range)]

        p = Pool(num_cores)
        p.map(get_extents, pool_args)
        p.terminate()

        num_frames = len(frame_range)
        this_azimuthal_extent_over_time = np.array(azimuthal_extent_over_time[:num_frames])
        this_radial_extent_over_time = np.array(radial_extent_over_time[:num_frames])
        this_radial_peak_over_time = np.array(radial_peak_over_time[:num_frames])

        #this_smoothed_differences_over_time = smooth(this_differences_over_time, 5)
        #this_growth_rates_over_time = np.diff(np.log(this_smoothed_differences_over_time)) / dt

        ##### Top Plot #####

        # Plot
        x = frame_range
        y1 = this_azimuthal_extent_over_time
        p1, = ax1.plot(x, y1, c = colors[i], linewidth = linewidth, zorder = 99 - i)

        # Axes
        if i == 3:
            ax1.set_xlim(x[0], x[-1])

        angles = np.linspace(0, 360, 7)
        ax1.set_ylim(0, 360)
        ax1.set_yticks(angles)

        # Annotate
        #ax1.set_xlabel("", fontsize = fontsize)
        ax1.set_ylabel("Azimuthal Extent (degrees)", fontsize = fontsize)

        alpha_coefficent = "3"
        if scale_height == 0.08:
            alpha_coefficent = "1.5"
        elif scale_height == 0.04:
            alpha_coefficent = "6"

        title1 = r"$h = %.2f$    $\alpha = %s \times 10^{%d}$" % (scale_height, alpha_coefficent, int(round(np.log(viscosity) / np.log(10), 0)) + 2)
        #title1 = r"$A = %.2f$" % (accretion)
        ax1.set_title("%s" % (title1), y = 1.035, fontsize = fontsize + 2)

        ##### Middle Plot #####
        y2 = this_radial_extent_over_time
        p2, = ax2.plot(x, y2, c = colors[i], linewidth = linewidth, zorder = 99 - i)

        # Axes
        if i == 3:
            ax2.set_xlim(x[0], x[-1])

        if scale_height == 0.08:
            ax2.set_ylim(0, 0.75)
        else:
            ax2.set_ylim(0, 0.50)

        # Annotate
        ax2.set_ylabel(r"Radial Extent (planet radii)", fontsize = fontsize)

        ##### Bottom Plot #####

        # Plot
        y3 = this_radial_peak_over_time
        p3, = ax3.plot(x, y3, c = colors[i], linewidth = linewidth, alpha = alpha, zorder = 90 - i)

        # Axes
        if i == 3:
            ax3.set_xlim(x[0], x[-1])

        if scale_height == 0.08:
            ax3.set_ylim(1.0, 2.5)
        else:
            ax3.set_ylim(1.0, 1.6)

        # Annotate
        ax3.set_xlabel("Time (planet orbits)", fontsize = fontsize)
        ax3.set_ylabel(r"Radial Center (planet radii)", fontsize = fontsize)
        #ax3.legend(loc = "upper right", fontsize = fontsize - 5)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/radiiAndExtentsOverTime_choice%d.png" % (save_directory, args.choice)
    else:
        save_fn = "%s/v%04d_radiiAndExtentsOverTime_choice%d.png" % (save_directory, version, arg.choice)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
