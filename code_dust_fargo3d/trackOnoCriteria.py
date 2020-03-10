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

    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "criteria",
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
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 17,
                         help = 'fontsize of plot annotations (default: 17)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 15,
                         help = 'fontsize of plot annotations (default: 15)')
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

### Data ###

density_zero = fromfile("gasdens0.dat").reshape(num_rad, num_theta) / surface_density_zero
averagedDensity_zero = np.average(density_zero, axis = 1)

def get_criteria(args_here):
    # Unwrap Args
    i, frame = args_here

    # Get Data
    density = fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta) / surface_density_zero
    averagedDensity = np.average(density, axis = 1)
    diff_averagedDensity = averagedDensity - averagedDensity_zero

    start_i = np.searchsorted(rad, 1.1)
    end_i = np.searchsorted(rad, 2.1)

    zoom_diff_averagedDensity = diff_averagedDensity[start_i : end_i]
    amplitude_i = np.argmax(zoom_diff_averagedDensity)

    amplitude = np.max(zoom_diff_averagedDensity)
    half_amplitude = np.e**(-0.5) * amplitude

    # Find half-max (left and right)
    left_zoom = np.abs(zoom_diff_averagedDensity[:amplitude_i] - half_amplitude)
    right_zoom = np.abs(zoom_diff_averagedDensity[amplitude_i:] - half_amplitude)

    left_i = np.argmin(left_zoom); right_i = np.argmin(right_zoom)
    left_r = rad[start_i + left_i]; right_r = rad[start_i + amplitude_i + right_i]; amplitude_r = rad[start_i + amplitude_i]
    width = (right_r - left_r) / amplitude_r

    amplitude_over_time[i] = amplitude / averagedDensity_zero[start_i + amplitude_i]
    width_over_time[i] = width / 2.0

    print i, frame, amplitude_over_time[i], width_over_time[i], left_r, amplitude_r, right_r


## Use These Frames ##
rate = 1 # 5 works better, but is very slow
start = 50
max_frame = 100 #util.find_max_frame()
#frame_range = np.array(range(start, max_frame + 1, rate))

amplitude_over_time = mp_array("d", len(frame_range))
width_over_time = mp_array("d", len(frame_range))
#contrasts_over_time = mp_array("d", len(frame_range))

for i, frame in enumerate(frame_range):
    get_criteria((i, frame))

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
    fig, host = plot.subplots()
    #fig.subplots_adjust(right = 0.75)
    #fig.subplots_adjust(right = 0.65)

    par1 = host.twinx()
    #par2 = host.twinx()

    # Plot
    x = frame_range
    y1 = amplitude_over_time
    y2 = width_over_time

    #ref, = par2.plot([x[0], x[-1]], [1.6, 1.6], c = 'k', linewidth = linewidth - 1) # to compare to Lindblad resonances (which we showed was useless)

    p1, = host.plot(x, y1, c = 'k', linewidth = linewidth)
    p2, = par1.plot(x, y2, c = 'b', linewidth = linewidth)

    #p4, = par4.plot(x, y3a, c = 'r', linewidth = linewidth)

    #p3, = par2.plot(x, y3, c = 'g', linewidth = linewidth)

    # Axes
    host.set_ylim(0, 1.1 * max(y1))
    par1.set_ylim(0, 1.6 * max(y2))

    #min_mass = args.min_mass; max_mass = args.max_mass; delta_mass = args.delta_mass
    #mass_ticks = np.arange(min_mass, max_mass, delta_mass)

    def tick_function(masses):
        # For the secondary x-axis showing the planet mass over time
        tick_locations = np.zeros(len(masses))
        tick_labels = []

        for i, mass in enumerate(masses):
            total_mass_jupiter = total_mass / jupiter_mass # in Jupiter masses
            times_i = az.my_searchsorted(total_mass_jupiter, mass)

            print mass, times_i, len(times)

            tick_locations[i] = times[times_i]
            if delta_mass < 0.1:
                tick_labels.append("%.2f" % mass)
            else:
                tick_labels.append("%.1f" % mass)

        return tick_locations, tick_labels

    #tick_locations, tick_labels = tick_function(mass_ticks)

    #par3.set_xlim(host.get_xlim())
    #par3.set_xticks(tick_locations)
    #par3.set_xticklabels(tick_labels)

    host.set_xlabel("Time (planet orbits)", fontsize = fontsize)
    host.set_ylabel(r"$\Sigma$ $/$ $\Sigma_0$", fontsize = fontsize)
    par1.set_ylabel(r"$\Delta w$ $/$ $r$", fontsize = fontsize, rotation = 270, labelpad = 15)
    #par2.set_ylabel("Radial Center (planet radii)", fontsize = fontsize, rotation = 270, labelpad = 20)
    #par3.set_xlabel(r"$M_\mathrm{p}$ [$M_\mathrm{J}$]", fontsize = fontsize)
    #par4.set_ylabel("Contrast", fontsize = fontsize, rotation = 270, labelpad = 20)

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    title1 = r"$h = %.2f$     $\alpha = %s \times 10^{%d}$     $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    plot.title("%s" % (title1), y = 1.035, fontsize = fontsize + 1)

    # Annotate
    tkw = dict(size=4, width=1.5)
    host.tick_params(axis = 'y', colors = p1.get_color(), **tkw)
    par1.tick_params(axis = 'y', colors = p2.get_color(), **tkw)
    #par2.tick_params(axis = 'y', colors = p3.get_color(), **tkw)
    #par3.tick_params(axis = 'x', **tkw)
    #par4.tick_params(axis = 'y', colors = p4.get_color(), **tkw)
    host.tick_params(axis = 'x', **tkw)

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1].split("-")[0]

    if version is None:
        save_fn = "%s/%s_criteriaOverTime.png" % (save_directory, directory_name)
    else:
        save_fn = "%s/v%04d_%s_criteriaOverTime.png" % (save_directory, version, directory_name)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)