"""
plot 2-D gas and dust density maps side-by-side

usage: plotBothDensityMaps.py [-h] [-c NUM_CORES] [--dir SAVE_DIRECTORY]
                              [--hide] [-v VERSION] [--range R_LIM R_LIM]
                              [--shift] [--cmap CMAP] [--cmaxGas GAS_CMAX]
                              [--cmaxDust DUST_CMAX] [--fontsize FONTSIZE]
                              [--dpi DPI]
                              frames [frames ...]

positional arguments:
  frames                select single frame or range(start, end, rate). error
                        if nargs != 1 or 3

optional arguments:
  -h, --help            show this help message and exit
  -c NUM_CORES          number of cores (default: 1)
  --dir SAVE_DIRECTORY  save directory (default: bothDensityMaps)
  --hide                for single plot, do not display plot (default: display
                        plot)
  -v VERSION            version number (up to 4 digits) for this set of plot
                        parameters (default: None)
  --range R_LIM R_LIM   radial range in plot (default: [r_min, r_max])
  --shift               center frame on vortex peak or middle (default: do not
                        center)
  --cmap CMAP           color map (default: viridis)
  --cmaxGas GAS_CMAX    maximum density in colorbar (default: 10 for hcm+, 2.5
                        otherwise)
  --cmaxDust DUST_CMAX  maximum density in colorbar (default: 10 for hcm+, 2.5
                        otherwise)
  --fontsize FONTSIZE   fontsize of plot annotations (default: 16)
  --dpi DPI             dpi of plot annotations (default: 100)
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
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
surface_density_zero = fargo_par["Sigma0"]
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
taper = fargo_par["MassTaper"]

size = fargo_par["PSIZE"]

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
    fig = plot.figure(figsize = (1400 / dpi, 600 / dpi), dpi = dpi)

    # Data
    gas_density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    dust_density = (fromfile("gasddens%d.dat" % frame).reshape(num_rad, num_theta))
    if center:
        if taper < 10.1:
            shift_c = az.get_azimuthal_peak(dust_density, fargo_par)
        else:
            shift_c = az.get_azimuthal_center(dust_density, fargo_par)
        gas_density = np.roll(gas_density, shift_c)
        dust_density = np.roll(dust_density, shift_c)

    ############################ Gas Density ##################################
    plot.subplot(1, 2, 1)

    # Data
    normalized_density = gas_density / surface_density_zero

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi)
    result = plot.pcolormesh(x, y, np.transpose(normalized_density), cmap = cmap)

    fig.colorbar(result)
    result.set_clim(gas_clim[0], gas_clim[1])

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    title = readTitle()

    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)

    if title is None:
        plot.title("Gas Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    else:
        plot.title("Gas Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)

    ############################ Dust Density #################################
    plot.subplot(1, 2, 2)

    # Data
    normalized_density = dust_density / (surface_density_zero / 100.0)

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi)
    result = plot.pcolormesh(x, y, np.transpose(normalized_density), cmap = cmap)

    fig.colorbar(result)
    result.set_clim(dust_clim[0], dust_clim[1])

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    title = readTitle()

    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)

    if title is None:
        plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    else:
        plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(0, 360)

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)

    ###########################################################################

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/bothDensityMaps_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_bothDensityMaps_%04d.png" % (save_directory, version, frame)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

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



