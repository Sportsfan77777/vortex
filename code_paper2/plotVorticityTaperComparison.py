"""
plot 2-D density maps
python plotDensityMaps.py
python plotDensityMaps.py frame_number
python plotDensityMaps.py -1 <<<===== Plots a sample
python plotDensityMaps.py -m
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
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import gridspec

from pylab import rcParams
from pylab import fromfile

import util
import utilVorticity
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

def new_argument_parser(description = "Plot dust density maps for four grain sizes in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = 2,
                         help = 'select two frames to compare intensity (first is for T = 10, second is for T = 1000)')

    # Directory Selection
    parser.add_argument('--dir1', dest = "directory1", default = '../taper10',
                         help = 'select first directory to compare vorticity (first is for T = 10, second is for T = 10) (default: ../taper10)')
    parser.add_argument('--dir2', dest = "directory2", default = '../taper1000',
                         help = 'select second directory to compare vorticity (first is for T = 10, second is for T = 1000) (default: ../taper1000)')

    parser.add_argument('-g', dest = "grain_size", default = "cm",
                         help = 'grain size string (default: cm)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "vorticityComparison",
                         help = 'save directory (default: vorticityComparison)')

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

    # Plot Parameters (contours)
    parser.add_argument('--contour', dest = "use_contours", action = 'store_true', default = False,
                         help = 'use contours or not (default: do not use contours)')
    parser.add_argument('--low', dest = "low_contour", type = float, default = 1.0,
                         help = 'lowest contour (default: 1.0)')
    parser.add_argument('--high', dest = "high_contour", type = float, default = 3.6,
                         help = 'highest contour (default: 3.5)')
    parser.add_argument('--num_levels', dest = "num_levels", type = int, default = None,
                         help = 'number of contours (choose this or separation) (default: None)')
    parser.add_argument('--separation', dest = "separation", type = float, default = 0.2,
                         help = 'separation between contours (choose this or num_levels) (default: 0.2)')

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "magma",
                         help = 'color map (default: magma)')
    parser.add_argument('--crange', dest = "c_lim", type = float, nargs = 2, default = [-0.2, 0],
                         help = 'range in colorbar (default: [-0.2, 0])')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'labelsize of plot annotations (default: 16)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get ID%04d Parameters ###
default_directory = "../taper1000/%s-size/" % (args.grain_size)
fargo_par = util.get_pickled_parameters(directory = default_directory)

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

### Get Input Parameters ###

# Frames
frame_range = args.frames

# Directories
grain_size = args.grain_size

directory1 = "%s/%s-size" % (args.directory1, grain_size)
directory2 = "%s/%s-size" % (args.directory2, grain_size)
directories = [directory1, directory2]

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Quantity to Plot
rossby = args.rossby
residual = args.residual

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
if args.c_lim is None:
    clim = [-0.2, 0]
else:
    clim = [args.c_lim[0], args.c_lim[1]]

fontsize = args.fontsize
labelsize = args.labelsize
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

def generate_colors(n):
    c = ['b', 'g', 'springgreen']
    colors = []
    for i in range(n):
        colors.append(c[i % len(c)])
    return colors

##### PLOTTING #####

def add_to_plot(frame, fig, ax, frame_i):
    # Declare Subplot
    #ax = plot.subplot(1, num_frames, frame_i, sharex = prev_ax, sharey = prev_ax, aspect = "equal")

    # Taper
    if frame_i == 1:
        taper_time = 10
    else:
        taper_time = 1000

    # Change directories
    cwd = os.getcwd()
    os.chdir(directories[frame_i - 1])

    # Data
    vrad = (fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta)) # add a read_vrad to util.py!
    vtheta = (fromfile("gasvtheta%d.dat" % frame).reshape(num_rad, num_theta)) # add a read_vrad to util.py!

    vorticity = utilVorticity.velocity_curl(vrad, vtheta, rad, theta, rossby = rossby, residual = residual)

    # Shift
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero
    dust_density = (fromfile("gasddens%d.dat" % frame).reshape(num_rad, num_theta))
    if center:
        if taper_time < 10.1:
            shift_c = az.get_azimuthal_peak(dust_density, fargo_par)
        else:
            threshold = util.get_threshold(size)
            shift_c = az.get_azimuthal_center(dust_density, fargo_par, threshold = threshold * surface_density_zero / 100.0)
        vorticity = np.roll(vorticity, shift_c)
        density = np.roll(density, shift_c)

    ### Plot ###
    x = rad
    y = theta * (180.0 / np.pi) - 180.0
    result = ax.pcolormesh(x, y, np.transpose(vorticity), cmap = cmap)

    #if frame_i == 2:
    #    cbar = fig.colorbar(result)

    result.set_clim(clim[0], clim[1])

    if use_contours:
        levels = np.linspace(low_contour, high_contour, num_levels)
        colors = generate_colors(num_levels)
        plot.contour(x, y, np.transpose(density), levels = levels, origin = 'upper', linewidths = 1, colors = colors, alpha = 0.8)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(-180, 180)

    angles = np.linspace(-180, 180, 7)
    plot.yticks(angles)

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    #title = readTitle()

    unit = "r_\mathrm{p}"
    plot.xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
    if frame_i == 1:
        plot.ylabel(r"$\phi - \phi_\mathrm{center}$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)
    plot.text(x_mid, y_text * (plot.ylim()[-1] - plot.ylim()[0]) + plot.ylim()[0], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_J$" % (int(planet_mass))
    text_visc = r"$\nu = 10^{%d}$" % (int(np.log(viscosity) / np.log(10)))
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    if frame_i == 1:
        plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * (plot.ylim()[-1] - plot.ylim()[0]) + plot.ylim()[0], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    if frame_i == 2:
        plot.text(0.84 * x_range / 2.0 + x_mid, y_text * (plot.ylim()[-1] - plot.ylim()[0]) + plot.ylim()[0], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')

    # Label colorbar
    # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
    if frame_i < 5:
        # Only for last frame
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "6%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)

        if frame_i == 2:
            if rossby:
               cbar_name = r"$\mathrm{Rossby}$ $\mathrm{number}$"
            else:
               cbar_name = r"$\mathrm{Vorticity}$"
            cbar.set_label(cbar_name, fontsize = fontsize, rotation = 270, labelpad = 30)

        #if frame_i != num_frames:
        #    fig.delaxes(cax) # to balance out frames that don't have colorbar with the one that does

    # Set Aspect Ratio
    #unit_aspect_ratio = (360) / (x_range)
    #ax.set_aspect(1.12 / unit_aspect_ratio)

    # Return to previous directory
    os.chdir(cwd)

    return ax
    
def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (14, 6), dpi = dpi)
    gs = gridspec.GridSpec(1, 2)

    frame_str = ""
    for i, frame_i in enumerate(frame_range):
        ax = fig.add_subplot(gs[i])
        ax = add_to_plot(frame_i, fig, ax, i + 1)
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    #### Finish Plot ####
    #title = r"$\mathrm{Beam:\ }\ \ %.03f^{\prime\prime} \times \ \ %.03f^{\prime\prime}$" % (arc_beam, arc_beam)
    #title = r"$%.03f^{\prime\prime} \times \ \ %.03f^{\prime\prime}$" % (arc_beam, arc_beam)
    #fig.suptitle(title, y = 0.945, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)

    # Save and Close
    plot.tight_layout()

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/vorticityTaperComparison_%s.png" % (save_directory, frame_str)
    else:
        save_fn = "%s/v%04d_vorticityTaperComparison_%s.png" % (save_directory, version, frame_str)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)