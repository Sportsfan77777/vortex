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
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pylab import rcParams
from pylab import fromfile

import util
import azimuthal as az
#from readTitle import readTitle

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
    parser.add_argument('frames', type = int, nargs = 2,
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "twoDustDensityMaps%d",
                         help = 'save directory (default: twoDustDensityMaps%d)')
    parser.add_argument('-n', dest = "dust_number", type = int, default = 1,
                         help = 'number (1, 2, or 3) corresponding to different dust sizes (default: 1)')
    parser.add_argument('--mpi', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpio output files (default: use dat)')
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

    # Plot Parameters (quiver)
    parser.add_argument('--quiver', dest = "quiver", action = 'store_true', default = False,
                         help = 'use velocity quivers or not (default: do not use quivers)')
    parser.add_argument('--start', dest = "start_quiver", type = float, default = None,
                         help = 'start of quiver region in radius (default: r_lim[0])')
    parser.add_argument('--end', dest = "end_quiver", type = float, default = None,
                         help = 'end of quiver region in radius (default: r_lim[1])')
    parser.add_argument('--rate_x', dest = "quiver_rate_x", type = int, default = 6,
                         help = 'sub_sample in radius (default: 6)')
    parser.add_argument('--rate_y', dest = "quiver_rate_y", type = int, default = 100,
                         help = 'sub_sample in angle (default: 24)')
    parser.add_argument('--scale', dest = "quiver_scale", type = float, default = 0.1,
                         help = 'bigger scale means smaller arrow (default: 0.1)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: inferno)')
    parser.add_argument('--cmax', dest = "cmax", type = float, default = 2,
                         help = 'maximum density in colorbar (default: 2)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 17,
                         help = 'fontsize of plot annotations (default: 17)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
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

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

scale_height = p.aspectratio
viscosity = p.nu

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
frame_range = args.frames

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory % args.dust_number
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

dust_number = args.dust_number
merge = args.merge
mpi = args.mpi

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

# Plot Parameters (quiver)
quiver = args.quiver
start_quiver = args.start_quiver
end_quiver = args.end_quiver
rate_x = args.quiver_rate_x
rate_y = args.quiver_rate_y
scale = args.quiver_scale
if start_quiver is None:
   start_quiver = x_min
if end_quiver is None:
   end_quiver = x_max

# Plot Parameters (constant)
cmap = args.cmap
clim = [0, args.cmax]

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

### Helper Functions ###

def find_min(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 0.9) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 1.25) # look for max radial density before r = 2.3
    min_rad_outer_index = np.argmin(averagedDensity[outer_disk_start : outer_disk_end])

    min_index = outer_disk_start + min_rad_outer_index
    min_rad = rad[min_index]
    min_density = averagedDensity[min_index]

    return min_density

def get_gap_depth(density):
    # Get Data
    averagedDensity = np.average(density, axis = 1)
    normalized_density = averagedDensity / surface_density_zero
    
    # Get Minima
    min_density = find_min(normalized_density)
    gap_depth =  1.0 / min_density

    return gap_depth

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

def make_plot(frames, show = False):
    # Set up figure
    fig = plot.figure(figsize = (8, 4), dpi = dpi)

    # Indvidual Subplots
    def add_to_plot(i):
        # Identify Subplot
        frame = frames[i]; number = i + 1
        ax = plot.subplot(1, 2, number)

        # Data
        field = "dens"
        if mpi:
          density = Fields("./", 'dust1', frame).get_field(field).reshape(num_z, num_rad, num_theta)
          gas_density = Fields("./", 'gas', frame).get_field(field).reshape(num_z, num_rad, num_theta)
        else:
          density = fromfile("dust1dens%d.dat" % frame).reshape(num_z, num_rad, num_theta)
          gas_density = fromfile("gasdens%d.dat" % frame).reshape(num_z, num_rad, num_theta)

        dz = z_angles[1] - z_angles[0]
        surface_density = np.sum(density[:, :, :], axis = 0) * dz
        gas_surface_density = np.sum(gas_density[:, :, :], axis = 0) * dz

        normalized_density = surface_density / dust_surface_density_zero # / np.sqrt(2.0 * np.pi) / scale_height_function[:, None]
        normalized_gas_density = gas_surface_density / surface_density_zero

        if center:
            normalized_density, shift_c = shift_density(normalized_density, fargo_par, reference_density = normalized_gas_density)
            normalized_gas_density, shift_c = shift_density(normalized_gas_density, fargo_par, reference_density = normalized_gas_density)

        ### Plot ###
        x = rad
        y = theta * (180.0 / np.pi)
        result = ax.pcolormesh(x, y, np.transpose(normalized_density), cmap = cmap)
        result.set_clim(clim[0], clim[1])

        # Contours
        if use_contours:
            levels = np.linspace(low_contour, high_contour, num_levels)
            colors = generate_colors(num_levels)
            plot.contour(x, y, np.transpose(normalized_gas_density), levels = levels, origin = 'upper', linewidths = 1, colors = colors)

        if quiver:
            # Velocity
            radial_velocity = np.array(fromfile("gasvy%d.dat" % frame).reshape(num_rad, num_theta)) # Radial
            azimuthal_velocity = np.array(fromfile("gasvx%d.dat" % frame).reshape(num_rad, num_theta)) # Azimuthal
            keplerian_velocity = rad * (np.power(rad, -1.5) - 1)
            azimuthal_velocity -= keplerian_velocity[:, None]

            if center:
                radial_velocity = np.roll(radial_velocity, shift_c, axis = -1)
                azimuthal_velocity = np.roll(azimuthal_velocity, shift_c, axis = -1)

            # Sub-sample the grid
            start_i = np.searchsorted(rad, start_quiver)
            end_i = np.searchsorted(rad, end_quiver)

            x_q = x[start_i:end_i]
            y_q = y[:]
            u = np.transpose(radial_velocity)[:, start_i:end_i]
            v = np.transpose(azimuthal_velocity)[:, start_i:end_i]

            plot.quiver(x_q[::rate_x], y_q[::rate_y], u[::rate_y,::rate_x], v[::rate_y,::rate_x], scale = scale, color = "cyan")

        # Axes
        plot.xlim(x_min, x_max)
        plot.ylim(0, 360)

        angles = np.linspace(0, 360, 7)
        plot.yticks(angles)

        # Annotate Axes
        orbit = (dt / (2 * np.pi)) * frame

        if orbit >= taper_time:
            current_mass = planet_mass
        else:
            current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

        current_mass += accreted_mass[frame]
        current_gap_depth = get_gap_depth(density)

        unit = "r_\mathrm{p}"
        plot.xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
        if number == 1:
           plot.ylabel(r"$\phi$ [degrees]", fontsize = fontsize)

        x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
        y_text = 1.14

        title = r"$t = %d$ [$m_\mathrm{p}=%.2f$ $M_\mathrm{J}$]" % (orbit, current_mass)
        #title = r"$t = %d$ [$\delta_\mathrm{gap}=%.1f$]" % (orbit, current_gap_depth)
        plot.title("%s" % (title), y = 1.035, fontsize = fontsize + 1)

        # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "6%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)
        cbar.set_label(r"Dust Surface Density  $\Sigma$ $/$ $\Sigma_0$", fontsize = fontsize, rotation = 270, labelpad = 25)

        if number != len(frames):
            fig.delaxes(cax) # to balance out frames that don't have colorbar with the one that does

    # Make each sub-plot
    for i, _ in enumerate(frames):
        add_to_plot(i)

    # Title
    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"
    title = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    #title = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $M_\mathrm{p} = %.2f$ $M_\mathrm{J}$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, planet_mass)
    plot.suptitle("%s" % (title), y = 1.10, fontsize = fontsize + 2, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1]

    if version is None:
        save_fn = "%s/%s_dustDensityMap_%04d-%04d.png" % (save_directory, directory_name, frames[0], frames[1])
    else:
        save_fn = "%s/v%04d_%s_dustDensityMap_%04d-%04d.png" % (save_directory, directory_name, version, frames[0], frames[1])
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi, pad_inches = 0.15)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####
make_plot(frame_range, show = show)

