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
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

def new_argument_parser(description = "Plot dust density maps for four grain sizes in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = 4,
                         help = 'select four frames to display the cm-size dust density maps')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "dustDensityEvolutionTwoByTwoMaps",
                         help = 'save directory (default: dustDensityEvolutionTwoByTwoMaps)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--box', dest = "box", type = float, default = 2.5,
                         help = 'width of box (in r_p) (default: 2.5)')
    parser.add_argument('--shift', dest = "center", action = 'store_false', default = True,
                         help = 'center frame on vortex peak or middle (default: center)')

    parser.add_argument('--cbar', dest = "colorbar", action = 'store_false', default = True,
                         help = 'include colorbar (default: no colorbar)')

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: magma)')
    parser.add_argument('--cmax', dest = "cmax", type = int, default = 300,
                         help = 'maximum density in colorbar (default: 300)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 18,
                         help = 'fontsize of plot annotations (default: 18)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 15,
                         help = 'fontsize of plot annotations (default: 15)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get ID%04d Parameters ###
fargo_par = util.get_pickled_parameters(directory = "../cm-size")

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

### Get Input Parameters ###

# Frames
frame_range = args.frames

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

version = args.version
box = args.box
center = args.center
colorbar = args.colorbar

# Plot Parameters (constant)
cmap = args.cmap
cmax = args.cmax
if cmax is not None:
    clim = [0, args.cmax]

fontsize = args.fontsize
labelsize = args.labelsize
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

def add_to_plot(frame, fig, ax, num_sizes, frame_i):
    # Convert size to number
    size_name = "cm"
    size = util.get_size(size_name)

    # Declare Subplot
    #ax = plot.subplot(1, num_frames, frame_i, sharex = prev_ax, sharey = prev_ax, aspect = "equal")

    # Data
    this_fargo_par = fargo_par.copy()
    this_fargo_par["PSIZE"] = size

    density = util.read_data(frame, 'dust', this_fargo_par, directory = "../%s-size" % size_name)

    if center:
        if taper_time < 10.1:
            shift = az.get_azimuthal_peak(density, fargo_par)
        else:
            threshold = util.get_threshold(size)
            shift = az.get_azimuthal_center(density, fargo_par, threshold = threshold * surface_density_zero)

        # Shift! (Handle gas case separately)
        density = np.roll(density, shift)
    normalized_density = density / surface_density_zero

    # Convert gas density to cartesian
    _, _, xs_grid, ys_grid, normalized_density = sq.polar_to_cartesian(normalized_density, rad, theta)

    ### Plot ###
    if size_name == "um":
        colormap = "viridis"
    else:
        colormap = cmap
    result = plot.pcolormesh(xs_grid, ys_grid, np.transpose(normalized_density), cmap = colormap)

    result.set_clim(clim[0], clim[1])

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    ax.add_artist(circle)

    # Add planet orbit
    planet_orbit = plot.Circle((0, 0), 1, color = "white", fill = False, alpha = 0.8, linestyle = "dashed", zorder = 50)
    ax.add_artist(planet_orbit)

    # Locate Planet
    if shift < -len(theta):
        shift += len(theta)
    planet_theta = theta[shift]
    planet_theta += (np.pi / 2.0) # Note: the conversion from polar to cartesian rotates everything forward by 90 degrees
    planet_theta = planet_theta % (2 * np.pi) # Keep 0 < theta < 2 * np.pi

    planet_x = np.cos(planet_theta)
    planet_y = np.sin(planet_theta)

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
    if frame_i > 2:
        ax.set_xlabel(r"$x$ [$r_p$]", fontsize = fontsize)
    if frame_i % 2 == 1:
        ax.set_ylabel(r"$y$ [$r_p$]", fontsize = fontsize)

    # Axes
    box_size = 2.5
    ax.set_xlim(-box_size, box_size)
    ax.set_ylim(-box_size, box_size)
    ax.set_aspect('equal')

    # Title
    title = r"$t$ $=$ $%.1f$  [$m_p(t)$ $=$ $%.2f$ $M_J$]" % (orbit, current_mass)
    plot.title("%s" % (title), y = 1.02, fontsize = fontsize + 1)

    # Title
    left_x = -0.8 * box_size; line_y = 1.15 * box_size; linebreak = 0.2 * box_size
    right_x = 1.3 * box_size
    if frame_i == 1:
        line1 = r'$M_p = %d$ $M_J$' % planet_mass
        line2 = r'$\nu = 10^{%d}$' % round(np.log(viscosity) / np.log(10), 0)
        plot.text(left_x, line_y + linebreak, line1, horizontalalignment = 'left', fontsize = fontsize + 2)
        plot.text(left_x, line_y, line2, horizontalalignment = 'left', fontsize = fontsize + 2)
    elif frame_i == 2 :
        line3 = r'$T_\mathrm{growth} = %d$ $\rm{orbits}$' % taper_time
        plot.text(right_x, line_y + 0.5 * linebreak, line3, horizontalalignment = 'right', fontsize = fontsize + 2)

    if frame_i <= 2:
        # Remove unless 1st frame
        ax.set_xticklabels([])

    # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
    if colorbar:
        # Only for last frame
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "8%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)
        if frame_i % 2 == 0:
            cbar.set_label(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{dust}$", fontsize = fontsize, rotation = 270, labelpad = 25)

        #if frame_i != num_frames:
        #    fig.delaxes(cax) # to balance out frames that don't have colorbar with the one that does

    return ax
    
def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (12, 10), dpi = dpi)
    gs = gridspec.GridSpec(2, 2)

    frame_str = ""
    for i, frame_i in enumerate(frame_range):
        ax = fig.add_subplot(gs[i])
        ax = add_to_plot(frame_i, fig, ax, len(frame_range), i + 1)
        frame_str += "%04d_" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    #### Finish Plot ####

    size = 1.0 # cm-size only (in the future, make this a parameter?)
    size_label = util.get_size_label(size)
    stokes_number = util.get_stokes_number(size)

    title = r"$\mathrm{1\ cm-size}$ $\mathrm{(St}_\mathrm{0}$ $=$ $%.03f \mathrm{)}$" % (stokes_number)
    fig.suptitle(title, y = 1.02, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)

    # Save and Close
    plot.tight_layout()

    # Save, Show,  and Close
    if version is None:
        save_fn = "%s/densityMapEvolution_%s.png" % (save_directory, frame_str)
    else:
        save_fn = "%s/v%04d_densityMapEvolution_%s.png" % (save_directory, version, frame_str)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)