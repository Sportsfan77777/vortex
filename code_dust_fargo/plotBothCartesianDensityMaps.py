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

def new_argument_parser(description = "Plot convolved intensity maps."):
    parser = argparse.ArgumentParser()

    # Directory Selection
    parser.add_argument('--dir1', dest = "directory1", default = None,
                         help = 'data directory for gas (default: None)')
    parser.add_argument('--dir2', dest = "directory2", default = None,
                         help = 'data directory for convolved cartesian intensity (default: None)')

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "sidebysideGasAndIntensity",
                         help = 'save directory (default: cartesianIntensityMaps)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--r_range', dest = "r_lim", type = int, nargs = 2, default = None,
                         help = 'id number for this set of plot parameters (default: [r_min, r_max])')
    parser.add_argument('-n', dest = "normalize", action = 'store_false', default = True,
                         help = 'normalize by max (default: normalize)')

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: inferno)')
    parser.add_argument('--cmax', dest = "cmax", type = int, default = None,
                         help = 'maximum density in colorbar (default: 2.5)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

dir1 = args.directory1
dir2 = args.directory2

### Get ID%04d Parameters ###
fn = "%s/id%04d_par.p" % (dir2, args.id_number)
fargo_par = pickle.load(open(fn, "rb"))

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

beam_size = fargo_par["Beam"]
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

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
normalize = args.normalize

# Plot Parameters (constant)
cmap = args.cmap
cmax = args.cmax
if cmax is not None:
    clim = [0, args.cmax]
elif normalize:
    cmax = 1
    clim = [0, 1]

fontsize = args.fontsize
dpi = args.dpi

###############################################################################

##### PLOTTING #####

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (14, 6), dpi = dpi)

    ############################## Left Panel #################################

    plot.subplot(1, 2, 1, aspect = 'equal')

    cwd = os.getcwd()
    os.chdir(dir1)

    # Data
    gas_fargo_par = util.get_pickled_parameters() ## shorten name?
    ######## Need to extract parameters, and add 'rad' and 'theta' ########
    gas_rad = np.linspace(gas_fargo_par['Rmin'], gas_fargo_par['Rmax'], gas_fargo_par['Nrad'])
    gas_theta = np.linspace(0, 2 * np.pi, gas_fargo_par['Nsec'])
    gas_fargo_par['rad'] = gas_rad; gas_fargo_par['theta'] = gas_theta
    gas_surface_density_zero = gas_fargo_par['Sigma0']

    gas_density = util.read_data(frame, 'gas', gas_fargo_par, id_number = id_number) / gas_surface_density_zero
    dust_density = util.read_data(frame, 'dust', gas_fargo_par, id_number = id_number)

    # Shift gas density with center of dust density
    shift = az.get_azimuthal_center(dust_density, gas_fargo_par, threshold = 10.0 * dust_surface_density_zero)
    gas_density = np.roll(gas_density, shift)

    # Locate Planet
    if shift < -len(gas_theta):
        shift += len(gas_theta)
    planet_theta = gas_theta[shift]
    planet_theta += (np.pi / 2.0) # Note: the conversion from polar to cartesian rotates everything forward by 90 degrees
    planet_theta = planet_theta % (2 * np.pi) # Keep 0 < theta < 2 * np.pi

    planet_x = np.cos(planet_theta)
    planet_y = np.sin(planet_theta)

    # Convert gas density to cartesian
    _, _, xs_grid, ys_grid, gas_density = sq.polar_to_cartesian(gas_density, gas_rad, gas_theta)

    ### Plot ###
    result = plot.pcolormesh(xs_grid, ys_grid, np.transpose(gas_density), cmap = "viridis")
    cbar = fig.colorbar(result)

    result.set_clim(0, 1.5)

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    fig.gca().add_artist(circle)

    # Add planet orbit
    planet_orbit = plot.Circle((0, 0), 1, color = "white", fill = False, alpha = 0.8, linestyle = "dashed", zorder = 50)
    fig.gca().add_artist(planet_orbit)

    # Label star and planet
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

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

    os.chdir(cwd)
    os.chdir(dir2)

    # Data
    intensity_cart = util.read_data(frame, 'cartesian_intensity', fargo_par, id_number = id_number)
    _, _, xs_grid, ys_grid = sq.get_cartesian_grid(rad)

    # Normalize
    if normalize:
        intensity_cart /= np.max(intensity_cart)

    ### Plot ###
    result = plot.pcolormesh(xs_grid, ys_grid, np.transpose(intensity_cart), cmap = cmap)
    cbar = fig.colorbar(result)

    if cmax is not None:
        result.set_clim(clim[0], clim[1])

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    fig.gca().add_artist(circle)

    # Add beam size
    beam = plot.Circle((-2, -2), beam_size / 2, color = "white")
    fig.gca().add_artist(beam)

    # Add planet orbit
    planet_orbit = plot.Circle((0, 0), 1, color = "white", fill = False, alpha = 0.8, linestyle = "dashed", zorder = 50)
    fig.gca().add_artist(planet_orbit)

    # Label star and planet
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    planet_size = current_mass / planet_mass
    plot.scatter(0, 0, c = "white", s = 300, marker = "*", zorder = 100) # star
    plot.scatter(planet_x, planet_y, c = "white", s = int(70 * planet_size), marker = "D", zorder = 100) # planet

    # Annotate Axes
    plot.xlabel("x", fontsize = fontsize)
    plot.ylabel("y", fontsize = fontsize)
    plot.title("Intensity Map (t = %.1f)" % (orbit), fontsize = fontsize + 1)

    # Axes
    box_size = 2.5
    plot.xlim(-box_size, box_size)
    plot.ylim(-box_size, box_size)

    ################################# End #####################################

    os.chdir(cwd)

    # Save, Show,  and Close
    if version is None:
        save_fn = "%s/id%04d_intensityCartGrid_%04d.png" % (save_directory, id_number, frame)
    else:
        save_fn = "%s/v%04d_id%04d_intensityCartGrid_%04d.png" % (save_directory, version, id_number, frame)
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



