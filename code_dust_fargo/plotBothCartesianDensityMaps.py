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
import azimuthal as az
import square as sq

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot dust density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "bothDensityMaps",
                         help = 'save directory (default: bothDensityMaps)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--shift', dest = "center", action = 'store_true', default = False,
                         help = 'center frame on vortex peak or middle (default: do not center)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "viridis",
                         help = 'color map (default: viridis)')
    parser.add_argument('--cmaxGas', dest = "gas_cmax", type = float, default = 2,
                         help = 'maximum density in colorbar (default: 2)')
    parser.add_argument('--cmaxDust', dest = "dust_cmax", type = float, default = None,
                         help = 'maximum density in colorbar (default: 10 for hcm+, 2.5 otherwise)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
fargo_par = util.get_pickled_parameters()

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
dust_surface_density_zero = surface_density_zero / 100
disk_mass = 2 * np.pi * dust_surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
center = args.center

# Plot Parameters (constant)
cmap = args.cmap
gas_cmax = args.gas_cmax
dust_cmax = args.dust_cmax
if dust_cmax is None:
    if size > 0.2:
        dust_cmax = 10
    else:
        dust_cmax = 2.5
gas_clim = [0, gas_cmax]
dust_clim = [0, dust_cmax]

fontsize = args.fontsize
dpi = args.dpi

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

##### PLOTTING #####

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (14, 6), dpi = dpi)

    # Data
    gas_density = util.read_data(frame, 'gas', fargo_par) / gas_surface_density_zero
    dust_density = util.read_data(frame, 'dust', fargo_par) / dust_surface_density_zero

    # Shift gas density with center of dust density
    threshold = util.get_threshold(size)
    shift = az.get_azimuthal_center(dust_density, fargo_par, threshold = threshold)
    gas_density = np.roll(gas_density, shift)
    dust_density = np.roll(dust_density, shift)

    # Locate Planet
    if shift < -len(theta):
        shift += len(theta)
    planet_theta = theta[shift]
    planet_theta += (np.pi / 2.0) # Note: the conversion from polar to cartesian rotates everything forward by 90 degrees
    planet_theta = planet_theta % (2 * np.pi) # Keep 0 < theta < 2 * np.pi

    planet_x = np.cos(planet_theta)
    planet_y = np.sin(planet_theta)

    # Convert density to cartesian
    _, _, xs_grid, ys_grid, gas_density = sq.polar_to_cartesian(gas_density, rad, theta)
    _, _, _, _, dust_density = sq.polar_to_cartesian(dust_density, rad, theta)

    ############################## Left Panel #################################

    plot.subplot(1, 2, 1, aspect = 'equal')

    ### Plot ###
    result = plot.pcolormesh(xs_grid, ys_grid, np.transpose(gas_density), cmap = "viridis")
    cbar = fig.colorbar(result)

    result.set_clim(dust_clim[0], dust_clim[1])

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    fig.gca().add_artist(circle)

    # Add planet orbit
    planet_orbit = plot.Circle((0, 0), 1, color = "white", fill = False, alpha = 0.8, linestyle = "dashed", zorder = 50)
    fig.gca().add_artist(planet_orbit)

    # Label star and planet
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    orbit = 550 + (time / (2 * np.pi)) * (frame - 550)
    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    # Use Fixed Mass
    current_mass = np.power(np.sin((np.pi / 2) * (550.0 / taper_time)), 2) * planet_mass

    planet_size = current_mass / planet_mass
    plot.scatter(0, 0, c = "white", s = 300, marker = "*", zorder = 100) # star
    plot.scatter(planet_x, planet_y, c = "white", s = int(70 * planet_size), marker = "D", zorder = 100) # planet

    # Annotate Axes
    plot.xlabel("x", fontsize = fontsize)
    plot.ylabel("y", fontsize = fontsize)
    plot.title("Gas Density Map (t = %.1f)" % (orbit), fontsize = fontsize + 1)

    # Axes
    box_size = 2.5
    plot.xlim(-box_size, box_size)
    plot.ylim(-box_size, box_size)

    ############################# Right Panel #################################

    plot.subplot(1, 2, 2, aspect = 'equal')

    ### Plot ###
    result = plot.pcolormesh(xs_grid, ys_grid, np.transpose(dust_density), cmap = cmap)
    cbar = fig.colorbar(result)

    result.set_clim(gas_clim[0], gas_clim[1])

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    fig.gca().add_artist(circle)

    # Add planet orbit
    planet_orbit = plot.Circle((0, 0), 1, color = "white", fill = False, alpha = 0.8, linestyle = "dashed", zorder = 50)
    fig.gca().add_artist(planet_orbit)

    # Label star and planet
    planet_size = current_mass / planet_mass
    plot.scatter(0, 0, c = "white", s = 300, marker = "*", zorder = 100) # star
    plot.scatter(planet_x, planet_y, c = "white", s = int(70 * planet_size), marker = "D", zorder = 100) # planet

    # Annotate Axes
    plot.xlabel("x", fontsize = fontsize)
    plot.ylabel("y", fontsize = fontsize)
    plot.title("Dust Density Map (t = %.1f)" % (orbit), fontsize = fontsize + 1)

    # Axes
    box_size = 2.5
    plot.xlim(-box_size, box_size)
    plot.ylim(-box_size, box_size)

    ################################# End #####################################

    # Save, Show,  and Close
    if version is None:
        save_fn = "%s/bothDensityMaps_%04d.png" % (save_directory, id_number, frame)
    else:
        save_fn = "%s/v%04d_bothDensityMaps_%04d.png" % (save_directory, version, id_number, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


###############################################################################

##### Make Plots! #####

# Iterate through frames

if len(frame_range) == 1:
    make_plot(frame_range[0], show = show)
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(make_plot, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            make_plot(frame)



