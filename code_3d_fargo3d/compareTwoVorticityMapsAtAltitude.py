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
import utilVorticity
import azimuthal as az
#from readTitle import readTitle

from advanced import Parameters
#from reader import Fields

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = 3,
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "twoAltitudeVorticityMaps",
                         help = 'save directory (default: fourVorticityMaps)')
    parser.add_argument('-m', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpio output files (default: use dat)')
    parser.add_argument('--merge', dest = "merge", type = int, default = 0,
                         help = 'number of cores needed to merge data outputs (default: 0)')

    # Quantity to plot
    parser.add_argument('--rossby', dest = "rossby", action = 'store_true', default = False,
                         help = 'plot rossby number instead of vorticity (default: plot vorticity)')
    parser.add_argument('--residual', dest = "residual", action = 'store_false', default = True,
                         help = 'use v_theta or v_theta - v_kep (default: use residual)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--shift', dest = "center", action = 'store_true', default = False,
                         help = 'center frame on vortex peak or middle (default: do not center)')

    parser.add_argument('-z', dest = "sliver", type = int, default = 0,
                         help = 'sliver above midplane to plot (default: 0)')

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
    parser.add_argument('--cmap', dest = "cmap", default = "magma",
                         help = 'color map (default: magma)')
    parser.add_argument('--crange', dest = "c_lim", type = float, nargs = 2, default = None,
                         help = 'range in colorbar (default: [-0.25, 0])')

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

num_rad = p.ny; num_theta = p.nx; num_z = p.nz
r_min = p.ymin; r_max = p.ymax
z_min = p.zmin; z_max = p.zmax

surface_density_zero = p.sigma0

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
frame = frame_range[0]
two_zs = frame_range[1:]

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

merge = args.merge
mpi = args.mpi

# Quantity to Plot
rossby = args.rossby
residual = args.residual

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)
z_angles = np.linspace(z_min, z_max, num_z)

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
if args.c_lim is None:
    clim = [-0.25, 0]
else:
    clim = [args.c_lim[0], args.c_lim[1]]

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

###############################################################################

def generate_colors(n):
    c = ['yellow', 'b', 'firebrick', 'w', 'green']
    colors = []
    for i in range(n):
        colors.append(c[i % len(c)])
    return colors

##### PLOTTING #####

def make_plot(frame, two_zs, show = False):
    # Set up figure
    fig = plot.figure(figsize = (8, 4), dpi = dpi)

    # Indvidual Subplots
    def add_to_plot(i):
        # Identify Subplot
        chosen_z = two_zs[i]; number = i + 1
        ax = plot.subplot(1, 2, number)

        # Data
        if mpi:
          density = Fields("./", 'gas', frame).get_field("dens").reshape(num_z, num_rad, num_theta)
          vrad = Fields("./", 'gas', frame).get_field("vy").reshape(num_z, num_rad, num_theta)
          vtheta = Fields("./", 'gas', frame).get_field("vx").reshape(num_z, num_rad, num_theta)
        else:
          density = fromfile("gasdens%d.dat" % frame).reshape(num_z, num_rad, num_theta)
          vrad = (fromfile("gasvy%d.dat" % frame).reshape(num_z, num_rad, num_theta)) # add a read_vrad to util.py!
          vtheta = (fromfile("gasvx%d.dat" % frame).reshape(num_z, num_rad, num_theta)) # add a read_vrad to util.py!
        midplane_density = density[num_z / 2 + chosen_z, :, :]
        midplane_vrad = vrad[num_z / 2 + chosen_z, :, :]
        midplane_vtheta = vtheta[num_z / 2 + chosen_z, :, :]

        dz = z_angles[1] - z_angles[0]
        surface_density = np.sum(density[:, :, :], axis = 0) * dz

        normalized_midplane_density = midplane_density / surface_density_zero # / np.sqrt(2.0 * np.pi) / scale_height_function[:, None]
        normalized_density = surface_density / surface_density_zero # / np.sqrt(2.0 * np.pi) / scale_height_function[:, None]

        vorticity = utilVorticity.velocity_curl(midplane_vrad, midplane_vtheta, rad, theta, rossby = rossby, residual = residual)

        ### Plot ###
        x = rad
        y = theta * (180.0 / np.pi)
        result = ax.pcolormesh(x, y, np.transpose(vorticity), cmap = cmap)
        result.set_clim(clim[0], clim[1])

        # Contours
        if use_contours:
            levels = np.linspace(low_contour, high_contour, num_levels)
            colors = generate_colors(num_levels)
            plot.contour(x, y, np.transpose(normalized_density), levels = levels, origin = 'upper', linewidths = 1, colors = colors)

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

        unit = "r_\mathrm{p}"
        plot.xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
        if number == 1:
           plot.ylabel(r"$\phi$ [degrees]", fontsize = fontsize)

        x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
        y_text = 1.14

        
        if i == 0:
          title = r"$\theta = \pi/2    t = %d$ [$m_\mathrm{p}=%.2f$ $M_\mathrm{J}$]" % (orbit, current_mass)
          plot.title("%s" % (title), y = 1.035, fontsize = fontsize + 1, loc = "left")
        else:
          this_z_angle = z_angles[num_z / 2 + z_level]
          this_z_angle = (this_z_angle - (np.pi / 2.0)) / scale_height
          title = r"$   \theta = \pi/2 + %.2f H$"
          plot.title("%s" % (title), y = 1.035, fontsize = fontsize + 1, loc = "right")

        # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "6%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)

        # Label colorbar
        if rossby:
           cbar_name = r"$\mathrm{Rossby}$ $\mathrm{number}$"
        else:
           cbar_name = r"$\mathrm{Vorticity}$"
        cbar.set_label(cbar_name, fontsize = fontsize, rotation = 270, labelpad = 25)

        if number != len(two_zs):
            fig.delaxes(cax) # to balance out frames that don't have colorbar with the one that does

    # Make each sub-plot
    for i, _ in enumerate(two_zs):
        add_to_plot(i)

    # Title
    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"
    #title = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    title = r"$\Sigma_0$ $/$ $\Sigma_\mathrm{base} = %.1f$    ($M_\mathrm{p} = %.2f$ $M_\mathrm{Jup}$)" % (surface_density_zero / surface_density_base, final_planet_mass)
    plot.suptitle("%s" % (title), y = 1.04, fontsize = fontsize + 2, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1]
    
    if version is None:
        save_fn = "%s/%s_vorticityMap_%04d-z%03d-z%03d.png" % (save_directory, directory_name, frame, two_zs[0], two_zs[1])
    else:
        save_fn = "%s/v%04d_%s_vorticityMap_%04d-z%03d-z%03d.png" % (save_directory, directory_name, version, frame, two_zs[0], two_zs[1])
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####
make_plot(frame, two_zs, show = show)

