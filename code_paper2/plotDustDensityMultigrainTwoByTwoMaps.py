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
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Size Selection
    parser.add_argument('--sizes', dest = "sizes", nargs = 4, default = ["um", "cm", "hcm", "mm"],
                         help = 'select 4 sizes (default: [cm, hcm, mm, um])')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "dustDensityMultigrainTwoByTwoMaps",
                         help = 'save directory (default: dustDensityMultigrainTwoByTwoMaps)')

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
    parser.add_argument('--cmax', dest = "cmax", type = int, default = 15,
                         help = 'maximum density in colorbar (default: 15), except for um (fixed: 1.5)')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
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
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Sizes
sizes = args.sizes

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

def add_to_plot(frame, fig, ax, size_name, num_sizes, frame_i):
    # Convert size to number
    size = util.get_size(size_name)

    # Declare Subplot
    #ax = plot.subplot(1, num_frames, frame_i, sharex = prev_ax, sharey = prev_ax, aspect = "equal")

    # Data
    this_fargo_par = fargo_par.copy()
    this_fargo_par["PSIZE"] = size

    if size_name == "um":
        # Gas case is separate!
        density = util.read_data(frame, 'dust', fargo_par, directory = "../cm-size")
        gas_density = util.read_data(frame, 'gas', fargo_par, directory = "../cm-size")
    else:
        density = util.read_data(frame, 'dust', this_fargo_par, directory = "../%s-size" % size_name)

    if center:
        if taper_time < 10.1:
            shift = az.get_azimuthal_peak(density, fargo_par)
        else:
            if size_name == "um":
                threshold = util.get_threshold(1.0) # get cm-size threshold instead
            else:
                threshold = util.get_threshold(size)
            shift = az.get_azimuthal_center(density, fargo_par, threshold = threshold * surface_density_zero)

        # Shift! (Handle gas case separately)
        if size_name == "um":
            density = np.roll(gas_density, shift) / 100.0
        else:
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

    if size_name == "um":
        result.set_clim(0, 1.5)
    else:
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
        ax.set_xlabel(r"$x$ [$r_\mathrm{p}$]", fontsize = fontsize)
    if frame_i % 2 == 1:
        ax.set_ylabel(r"$y$ [$r_\mathrm{p}$]", fontsize = fontsize)

    # Axes
    box_size = 2.5
    ax.set_xlim(-box_size, box_size)
    ax.set_ylim(-box_size, box_size)
    ax.set_aspect('equal')

    if frame_i > 1:
        ax.spines['bottom'].set_color('w'); ax.spines['top'].set_color('w'); ax.spines['left'].set_color('w'); ax.spines['right'].set_color('w')
        ax.tick_params(colors = 'white', labelcolor = 'black', width = 1, length = 5)

    # Label
    size_label = util.get_size_label(size)
    stokes_number = util.get_stokes_number(size)

    title = r"%s$\mathrm{-size}$" % size_label
    stokes = r"$\mathrm{St}_\mathrm{0}$ $=$ $%.03f$" % stokes_number
    if size_name == "um":
        title = r"$\mathrm{Gas\ Density}$"
        plot.text(-0.9 * box_size, 2, title, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    else:
        plot.text(-0.9 * box_size, 2, title, fontsize = fontsize, color = 'white', horizontalalignment = 'left', bbox=dict(facecolor = 'black', edgecolor = 'white', pad = 10.0))
        plot.text(0.9 * box_size, 2, stokes, fontsize = fontsize, color = 'white', horizontalalignment = 'right')
    #ax.set_title(title)
    frame_title = r"$t$ $=$ $%.1f$ [$m_p(t)$ $=$ $%.2f$ $M_J$]" % (orbit, current_mass)

    # Title
    left_x = -0.8 * box_size; line_y = 1.1 * box_size; linebreak = 0.2 * box_size
    right_x = 1.3 * box_size
    if frame_i == 1:
        line1 = r'$M_p = %d$ $M_J$' % planet_mass
        line2 = r'$\nu = 10^{%d}$' % round(np.log(viscosity) / np.log(10), 0)
        plot.text(left_x, line_y + linebreak, line1, horizontalalignment = 'left', fontsize = fontsize + 2)
        plot.text(left_x, line_y, line2, horizontalalignment = 'left', fontsize = fontsize + 2)

        if orbit < taper_time:
            # Add this white line to keep supertitle in the save frame
            plot.text(left_x, line_y + 2.0 * linebreak, line1, color = "white", horizontalalignment = 'left', fontsize = fontsize + 2)
    elif frame_i == 2 :
        line3 = r'$T_\mathrm{growth} = %d$ $\rm{orbits}$' % taper_time
        plot.text(right_x, line_y + 0.5 * linebreak, line3, horizontalalignment = 'right', fontsize = fontsize + 2)

    #if frame_i != 1:
    #    # Remove unless 1st frame
    #    ax.set_yticklabels([])

    # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
    if colorbar:
        # Only for last frame
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "8%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)
        if frame_i == 1:
            cbar.set_label(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{gas}$", fontsize = fontsize, rotation = 270, labelpad = 25)
        else:
            cbar.set_label(r"$\Sigma$ / $\Sigma_\mathrm{0,}$ $_\mathrm{dust}$", fontsize = fontsize, rotation = 270, labelpad = 25)

        #if frame_i != num_frames:
        #    fig.delaxes(cax) # to balance out frames that don't have colorbar with the one that does

    return ax
    
def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (12, 10), dpi = dpi)
    gs = gridspec.GridSpec(2, 2)

    size_str = ""
    for i, size_i in enumerate(sizes):
        ax = fig.add_subplot(gs[i])
        ax = add_to_plot(frame, fig, ax, size_i, len(frame_range), i + 1)
        size_str += "size_i_"
    size_str = size_str[:-1] # Trim last '_'

    #### Finish Plot ####

    # Title
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    current_mass = util.get_current_mass(orbit, taper_time, planet_mass = planet_mass)
    if orbit >= taper_time:
        frame_title = r"$t$ $=$ $%.1f$" % (orbit)
        fig.suptitle(frame_title, y = 1.0, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)
    else:
        frame_title = "\n" + r"$t$ $=$ $%.1f$" % (orbit) + "\n" + "[$m_p(t)$ $=$ $%.2f$ $M_J$]" % (current_mass)
        fig.suptitle(frame_title, y = 1.0, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 20.0), fontsize = fontsize + 4)

    # Save and Close
    plot.tight_layout()

    # Save, Show,  and Close
    png = "png"; pdf = "pdf"
    if version is None:
        save_fn = "%s/densityMap_%04d.%s" % (save_directory, frame, png)
        pdf_save_fn = "%s/densityMap_%04d.%s" % (save_directory, frame, pdf)
    else:
        save_fn = "%s/v%04d_densityMap_%04d.%s" % (save_directory, version, frame, png)
        pdf_save_fn = "%s/v%04d_densityMap_%04d.%s" % (save_directory, version, frame, pdf)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)
    #plot.savefig(pdf_save_fn, bbox_inches = 'tight', format = "pdf") # pdf will not work

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
