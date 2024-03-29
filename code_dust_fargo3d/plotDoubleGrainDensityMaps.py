"""
plot 2-D dust density maps

usage: plotDustDensityMaps.py [-h] [-c NUM_CORES] [--dir SAVE_DIRECTORY]
                              [--hide] [-v VERSION] [--range R_LIM R_LIM]
                              [--shift] [--cmap CMAP] [--cmax CMAX]
                              [--fontsize FONTSIZE] [--dpi DPI]
                              frames [frames ...]

positional arguments:
  frames                select single frame or range(start, end, rate). error
                        if nargs != 1 or 3

optional arguments:
  -h, --help            show this help message and exit
  -c NUM_CORES          number of cores (default: 1)
  --dir SAVE_DIRECTORY  save directory (default: dustDensityMaps)
  --hide                for single plot, do not display plot (default: display
                        plot)
  -v VERSION            version number (up to 4 digits) for this set of plot
                        parameters (default: None)
  --range R_LIM R_LIM   radial range in plot (default: [r_min, r_max])
  --shift               center frame on vortex peak or middle (default: do not
                        center)
  --cmap CMAP           color map (default: viridis)
  --cmax CMAX           maximum density in colorbar (default: 10 for hcm+, 2.5
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

from advanced import Parameters
from reader import Fields

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "double-grainDensityMaps",
                         help = 'save directory (default: multigrainDensityMaps)')
    parser.add_argument('--mpi', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpio output files (default: do not use dat)')
    parser.add_argument('--merge', dest = "merge", type = int, default = 0,
                         help = 'number of cores needed to merge data outputs (default: 0)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--shift', dest = "center", action = 'store_true', default = False,
                         help = 'center frame on vortex peak or middle (default: do not center)')

    # Plot Parameters (contours)
    parser.add_argument('--contour', dest = "use_contours", action = 'store_true', default = False,
                         help = 'use contours or not (default: do not use contours)')
    parser.add_argument('--low', dest = "low_contour", type = float, default = 1.1,
                         help = 'lowest contour (default: 1.1)')
    parser.add_argument('--high', dest = "high_contour", type = float, default = 3.5,
                         help = 'highest contour (default: 3.5)')
    parser.add_argument('--num_levels', dest = "num_levels", type = int, default = None,
                         help = 'number of contours (choose this or separation) (default: None)')
    parser.add_argument('--separation', dest = "separation", type = float, default = 0.1,
                         help = 'separation between contours (choose this or num_levels) (default: 0.1)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: inferno)')
    parser.add_argument('--cmax', dest = "cmax", type = float, nargs = 2, default = None,
                         help = 'dust maximum density in colorbar (default: None)')
    parser.add_argument('--cmaxGas', dest = "cmaxGas", type = float, default = 1.4,
                         help = 'gas maximum density in colorbar (default: 1.4)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 18,
                         help = 'fontsize of plot annotations (default: 18)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'labelsize of plot annotations (default: 16)')
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
dust_surface_density_zero = p.sigma0 * p.epsilon

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

viscosity = p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()

"""
num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
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

merge = args.merge
mpi = args.mpi

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
if args.r_lim is None:
    y_min = r_min; y_max = r_max
else:
    y_min = args.r_lim[0]; y_max = args.r_lim[1]
center = args.center

# Plot Parameters (contours)
use_contours = args.use_contours
low_contour = args.low_contour
high_contour = args.high_contour
num_levels = args.num_levels
if num_levels is None:
    separation = args.separation
    num_levels = int(round((high_contour - low_contour) / separation + 1, 0))

# Plot Parameters (constant)
cmap = args.cmap
cmax = args.cmax
cmaxGas = args.cmaxGas

fontsize = args.fontsize
labelsize = args.labelsize
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

# Planet File
# Data
data = np.loadtxt("planet0.dat")
times = data[:, 0]; base_mass = data[:, 7]
accreted_mass = data[:, 8] / jupiter_mass


"""
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

# Plot Parameters (contours)
use_contours = args.use_contours
low_contour = args.low_contour
high_contour = args.high_contour
num_levels = args.num_levels
if num_levels is None:
    separation = args.separation
    num_levels = int(round((high_contour - low_contour) / separation + 1, 0))

# Plot Parameters (constant)
cmap = args.cmap
clim = [0, args.cmax]

fontsize = args.fontsize
dpi = args.dpi
"""

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Helper Functions ###

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

###############################################################################

def generate_colors(n):
    c = ['yellow', 'b', 'firebrick', 'w', 'green']
    colors = []
    for i in range(n):
        colors.append(c[i % len(c)])
    return colors

##### PLOTTING #####

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 11), dpi = dpi)

    # Indvidual Subplots (1 Gas + 3 Dust)
    def add_to_plot(i, grain):
        # Identify Subplot
        number = i + 1
        ax = plot.subplot(3, 1, number)

        # Data
        if merge > 0:
            num_merged_cores = merge
            gas_density = util.read_merged_data(frame, num_merged_cores, num_rad, num_theta)
            if i > 0:
               density = util.read_merged_data(frame, num_merged_cores, num_rad, num_theta, fluid = 'dust%d' % i)
        elif mpi:
            field = "dens"
            gas_density = Fields("./", 'gas', frame).get_field(field).reshape(num_rad, num_theta)
            if i > 0:
               density = Fields("./", 'dust%d' % i, frame).get_field(field).reshape(num_rad, num_theta)
        else:
            gas_density = fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)
            if i > 0:
               density = fromfile("dust%ddens%d.dat" % (i, frame)).reshape(num_rad, num_theta)
        normalized_gas_density = gas_density / surface_density_zero
        if i > 0:
           normalized_density = density / dust_surface_density_zero

        if center:
            if i > 0:
                normalized_density, shift_c  = shift_density(normalized_density, fargo_par, reference_density = normalized_gas_density)
            normalized_gas_density, shift_c = shift_density(normalized_gas_density, fargo_par, reference_density = normalized_gas_density)

        ### Plot ###
        x = theta * (180.0 / np.pi)
        y = rad

        if i == 0:
           cmap = 'viridis'
           result = ax.pcolormesh(x, y, normalized_gas_density, cmap = cmap)
           cbar = fig.colorbar(result); result.set_clim(0, cmaxGas)
        else:
           cmap =  args.cmap
           result = ax.pcolormesh(x, y, normalized_density, cmap = cmap)
           cbar = fig.colorbar(result); result.set_clim(0, cmax[i - 1])

        if number == 1:
            cbar.set_label(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{gas}$", fontsize = fontsize, rotation = 270, labelpad = 25)
        elif number == 3:
            cbar.set_label(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{dust}$", fontsize = fontsize, rotation = 270, labelpad = 25)

        if use_contours and i > 0:
            levels = np.linspace(low_contour, high_contour, num_levels)
            colors = generate_colors(num_levels)
            plot.contour(x, y, normalized_gas_density, levels = levels, origin = 'upper', linewidths = 1, colors = colors)

        # Axes
        plot.xlim(0, 360)
        plot.ylim(y_min, y_max)
        
        angles = np.linspace(0, 360, 7)
        plot.xticks(angles)

        ax.spines['bottom'].set_color('w'); ax.spines['top'].set_color('w'); ax.spines['left'].set_color('w'); ax.spines['right'].set_color('w')
        ax.tick_params(colors = 'white', labelcolor = 'black', width = 1, length = 6)


        # Annotate Axes
        orbit = (dt / (2 * np.pi)) * frame

        if orbit >= taper_time:
            current_mass = planet_mass
        else:
            current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

        current_mass += accreted_mass[frame]

        if number == 3:
           plot.xlabel(r"$\phi$ $\mathrm{(degrees)}$", fontsize = fontsize)
        if number == 1 or number == 3:
           plot.ylabel(r"Radius [$r_\mathrm{p}$]", fontsize = fontsize)

        if number == 1:
           plot.title(r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{J}$]" % (orbit, current_mass), bbox=dict(facecolor = 'w', edgecolor = 'k', pad = 10.0), y = 1.08, fontsize = fontsize + 1)

        # Label
        left_x = plot.xlim()[0]; right_x = plot.xlim()[-1]; range_x = right_x - left_x; margin_x = 0.05 * range_x
        bottom_y = plot.ylim()[0]; top_y = plot.ylim()[-1]; range_y = top_y - bottom_y; margin_y = 0.15 * range_y

        if number == 1:
           # Gas
           title = r"$\mathrm{Gas\ Density}$"
           #plot.text(left_x + margin_x, top_y - margin_y, title, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
        else:
           # Dust
           this_size = util.get_size(grain)

           size_label = util.get_size_label(this_size)
           stokes_number = util.get_stokes_number(this_size)

           title = r"%s$\mathrm{-size}$" % size_label
           stokes = r"$\mathrm{St}_\mathrm{0}$ $=$ $%.03f$" % stokes_number
           
           #plot.text(left_x + margin_x, top_y - margin_y, title, fontsize = fontsize, color = 'white', horizontalalignment = 'left', bbox=dict(facecolor = 'black', edgecolor = 'white', pad = 10.0))
           #plot.text(right_x - margin_x, top_y - margin_y, stokes, fontsize = fontsize, color = 'white', horizontalalignment = 'right', bbox=dict(facecolor = 'black', edgecolor = 'white', pad = 10.0))

        # Text
        line_y = top_y + 0.31 * range_y; linebreak = 0.16 * range_y
        left_start_x = left_x - 3.0 * margin_x; right_end_x = right_x + 4.0 * margin_x;
        if number == 1:
           line1 = r'$M_p = %d$ $M_J$' % planet_mass
           line2 = r'$\nu = 10^{%d}$' % round(np.log(viscosity) / np.log(10), 0)
           plot.text(left_start_x, line_y + linebreak, line1, horizontalalignment = 'left', fontsize = fontsize + 1)
           plot.text(left_start_x, line_y, line2, horizontalalignment = 'left', fontsize = fontsize + 1)

           line3 = r'$T_\mathrm{growth} = %d$ $\rm{orbits}$' % taper_time
           line4 = r"$N_\mathrm{r} \times \ N_\mathrm{\phi} = %d \times \ %d$" % (num_rad, num_theta)

           plot.text(right_end_x, line_y + linebreak, line3, horizontalalignment = 'right', fontsize = fontsize + 1)
           plot.text(right_end_x, line_y, line4, horizontalalignment = 'right', fontsize = fontsize + 1)

    add_to_plot(0, "gas") # Actually Gas!
    add_to_plot(1, "cm")
    add_to_plot(2, "hcm")

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/multigrainDensityMaps_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_multigrainDensityMaps_%04d.png" % (save_directory, version, frame)
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

