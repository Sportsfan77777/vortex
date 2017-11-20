"""
convolves intensity map to make it look like an alma image

python convolveIntensityMap.py frame
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot convolved intensity maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "intensityMaps",
                         help = 'save directory (default: intensityMaps)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('-s', dest = "new_res", nargs = 2, type = int, default = [400, 400],
                         help = 're-sample resolution (default: [400, 400])')
    parser.add_argument('--r_range', dest = "r_lim", type = int, nargs = 2, default = None,
                         help = 'id number for this set of plot parameters (default: [r_min, r_max])')
    

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "viridis",
                         help = 'color map (default: viridis)')
    parser.add_argument('--cmax', dest = "cmax", type = int, default = 2.5,
                         help = 'maximum density in colorbar (default: 2.5)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get ID%04d Parameters ###
fn = "id%04d_par.p" % args.id_number
fargo_par = pickle.load(open(fn, "rb"))

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

beam = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

### Get Input Parameters ###

# Frames
if len(args.frames) == 1:
    frame_range = args.frames
elif len(args.frames) == 3:
    start = args.frames[0]; end = args.frames[1]; rate = args.frames[2]
    frame_range = range(start, end + 1, rate)
else:
    print "Error: Must supply 1 or 3 frame arguments\nWith one argument, plots single frame\nWith three arguments, plots range(start, end + 1, rate)"
    exit()

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show

new_num_rad = args.new_res[0]; new_num_theta = args.new_res[1]
rad = np.linspace(r_min, r_max, new_num_rad)
theta = np.linspace(0, 2 * np.pi, new_num_theta)

id_number = args.id_number
version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]

# Plot Parameters (constant)
cmap = args.cmap
clim = [0, args.cmax]

fontsize = args.fontsize
dpi = args.dpi

###############################################################################

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "inferno"
#clim = [0, 2]
clim_in = np.searchsorted(rad, 1.1)
clim_out = np.searchsorted(rad, 2.3)

clim_low = 10
clim_high = 96

fontsize = 14
my_dpi = 100

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    intensity = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number)

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi)
    result = ax.pcolormesh(x, y, np.transpose(intensity), cmap = cmap)

    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)
    plot.title("Intensity Map (t = %.1f)" % (orbit), fontsize = fontsize + 1)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)

    # Save, Show,  and Close
    if version is None:
        save_fn = "%s/id%04d_intensityMap%04d.png" % (save_directory, id_number, frame)
    else:
        save_fn = "%s/v%04d_id%04d_intensityMap%04d.png" % (save_directory, version, id_number, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


###############################################################################

def old_make_plot(frame, show = False):
    """
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    print frame

    ##### Get rid of choose_axis (just do one axis!) #####

    def choose_axis(i, axis):
        # Orbit Number
        time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
        orbit = int(round(time / (2 * np.pi), 0)) * i

        # Mass
        if orbit >= taper_time:
            current_mass = planet_mass / 0.001
        else:
            current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * (planet_mass / 0.001)

        # Set up figure
        fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
        ax = fig.add_subplot(111)

        # Data  <<<---- Read in data outside of this function.
        data = np.loadtxt("intensitymap.out")
        intensity = data[:, -1].reshape(num_rad, num_theta)
        intensity = clear_inner_disk(intensity)
        xs, ys, xs_grid, ys_grid, intensity_cart = polar_to_cartesian(intensity, rad, theta)

        # Convolve Data
        convolved_intensity = convolve_intensity(intensity_cart)
        contrast, maximum, opposite = record_contrast(convolved_intensity, xs, ys)

        if axis == "zoom":
            prefix = "zoom_"
            sq = 2.5
        else:
            prefix = ""
            sq = np.max(xs_grid)

        plot.xlim(-sq, sq)
        plot.ylim(-sq, sq)
        plot.axes().set_aspect('equal')

        ### Plot ###
        result = ax.pcolormesh(xs_grid, ys_grid, np.transpose(convolved_intensity), cmap = cmap)
        cbar = fig.colorbar(result)

        clim = [np.percentile(intensity[clim_in : clim_out, :], clim_low), np.percentile(intensity[clim_in : clim_out, :], clim_high)]
        result.set_clim(clim[0], clim[1])

        cbar.set_label(r"Intensity", fontsize = fontsize + 2, rotation = 270, labelpad = 20)

        # Get rid of interior
        circle = plot.Circle((0, 0), min(rad), color = "black")
        fig.gca().add_artist(circle)

        # Add beam size
        beam = plot.Circle((-2, -2), beam_size, color = "white")
        fig.gca().add_artist(beam)

        # Label star and planet
        planet_size = (current_mass / (planet_mass / 0.001))
        plot.scatter(0, 0, c = "white", s = 300, marker = "*", zorder = 100) # star
        plot.scatter(0, 1, c = "white", s = int(70 * planet_size), marker = "D") # planet

        # Add minor grid lines
        alpha = 0.2
        dashes = [1, 5]
        #plot.grid(b = True, which = 'major', color = "black", dashes = dashes, alpha = alpha)
        #plot.grid(b = True, which = 'minor', color = "black", dashes = dashes, alpha = alpha)
        plot.minorticks_on()

        # Annotate
        title1 = r"$m_p = %d$ $M_{Jup}$, $\nu_{disk} = 10^{%d}$, $T_{growth} = %d$ $\rm{orbits}$" % (int(planet_mass / 0.001), int(np.log(viscosity) / np.log(10)), taper_time)
        title2 = r"$t = %d$ $\rm{orbits}}$, $m_p(t) = %.2f$ $M_{Jup}$ | $r = %d$ $\rm{AU}$, $\lambda = %d$ $\rm{\mu m}$" % (orbit, current_mass, int(radius), wavelength)
        #plot.xlabel("x", fontsize = fontsize)
        #plot.ylabel("y", fontsize = fontsize)
        plot.title("%s" % (title2), y = 1.01, fontsize = fontsize)
        plot.text(0.0, 3.14, title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

        # Write contrast
        plot.text(-2, 2, "Contrast: %.2f = %.2e / %.2e" % (contrast, maximum, opposite), color = "white", fontsize = fontsize)

        # Save and Close
        plot.savefig("%sconvolvedIntensityMap%04d_r%d_at%d.png" % (prefix, i, int(radius), wavelength), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    #choose_axis(i, "normal")
    choose_axis(i, "zoom")
    """


##### Make Plots! #####

# Iterate through frames

if len(frame_range) == 1:
    make_plot(frame_range[0])
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(make_plot, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            make_plot(frame)



