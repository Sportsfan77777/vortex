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
    parser.add_argument('frames', type = int, nargs = 2,
                         help = 'select two frames to compare intensity (first is for T = 10, second is for T = 1000)')

    # Directory Selection
    parser.add_argument('--dir1', dest = "directory1", default = '../taper10',
                         help = 'select first directory to compare intensity (first is for T = 10, second is for T = 10) (default: ../taper10)')
    parser.add_argument('--dir2', dest = "directory2", default = '../taper750',
                         help = 'select second directory to compare intensity (first is for T = 10, second is for T = 1000) (default: ../taper1000)')

    parser.add_argument('-w', dest = "wavelength", type = float, default = 870,
                         help = 'wavelength (in um) (default: 870)')
    parser.add_argument('-b', dest = "beam_size", type = float, default = 10,
                         help = 'beam_size (in AU) (default: 10)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "intensityComparison",
                         help = 'save directory (default: intensityComparison)')

    # Old Format
    parser.add_argument('--old_res', dest = "old_res", type = int, nargs = 2, default = [1024, 2048],
                         help = 'select two frames to display the intensity maps')


    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--box', dest = "box", type = float, default = 2.5,
                         help = 'width of box (in r_p) (default: 2.5)')
    parser.add_argument('--arc', dest = "arc", action = 'store_false', default = True,
                         help = 'axes in arcseconds (default: yes, arcseconds!)')
    parser.add_argument('-n', dest = "normalize", action = 'store_false', default = True,
                         help = 'normalize by max (default: normalize)')

    parser.add_argument('--cbar', dest = "colorbar", action = 'store_false', default = True,
                         help = 'include colorbar (default: no colorbar)')

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: magma)')
    parser.add_argument('--cmax', dest = "cmax", type = int, default = None,
                         help = 'maximum density in colorbar (default: None, except 1 if normalized)')

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
default_directory = "taper750/synthetic/lambda%04d/beam%03d" % (args.wavelength, args.beam_size)

fn = "../%s/id%04d_par.p" % (default_directory, args.id_number)
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

planet_radius = fargo_par["Radius"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

### Get Input Parameters ###

# Frames
frame_range = args.frames

# Directories
directory1 = "%s/synthetic" % args.directory1
directory2 = "%s/synthetic" % args.directory2
directories = [directory1, directory2]

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Old Format
old_num_rad = args.old_res[0]
old_num_theta = args.old_res[1]

# Plot Parameters (variable)
show = args.show

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

box = args.box
arc = args.arc
normalize = args.normalize
colorbar = args.colorbar

# Plot Parameters (constant)
cmap = args.cmap
cmax = args.cmax
if cmax is not None:
    clim = [0, args.cmax]
elif normalize:
    cmax = 1
    clim = [0, 1]

fontsize = args.fontsize
labelsize = args.labelsize
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

#directory_base = "../taper%d/synthetic"
#directories = [directory_base % 10, directory_base % 1000]

def add_to_plot(frame, fig, ax, num_sizes, frame_i):
    # Declare Subplot
    #ax = plot.subplot(1, num_frames, frame_i, sharex = prev_ax, sharey = prev_ax, aspect = "equal")

    # Taper
    if frame_i == 1:
        taper_time = 10
    else:
        taper_time = 750

    # Change directories
    cwd = os.getcwd()
    os.chdir(directories[frame_i - 1])

    # Data
    intensity_cart = util.read_data(frame, 'cartesian_intensity', fargo_par, id_number = id_number, directory = "lambda%04d/beam%03d" % (args.wavelength, args.beam_size))
    xs, ys, xs_grid, ys_grid = sq.get_cartesian_grid(rad)

    # Get Shift
    dust_fargo_par = util.get_pickled_parameters(directory = "../cm-size") ## shorten name?
    ######## Need to extract parameters, and add 'rad' and 'theta' ########
    dust_rad = np.linspace(dust_fargo_par['Rmin'], dust_fargo_par['Rmax'], dust_fargo_par['Nrad'])
    dust_theta = np.linspace(0, 2 * np.pi, dust_fargo_par['Nsec'])
    dust_fargo_par['rad'] = dust_rad; dust_fargo_par['theta'] = dust_theta
    gas_surface_density_zero = dust_fargo_par['Sigma0']

    dust_density = util.read_data(frame, 'dust', dust_fargo_par, id_number = id_number, directory = "../cm-size")

    # Shift gas density with center of dust density
    shift = az.get_azimuthal_center(dust_density, dust_fargo_par, threshold = 10.0 * gas_surface_density_zero / 100.0)

    # Normalize
    if normalize:
        intensity_cart /= np.max(intensity_cart)

    # Arcseconds or Planet Radii
    if arc:
        arc_weight = planet_radius / distance # related to parallax
    else:
        arc_weight = 1

    ## Plot! ##
    result = plot.pcolormesh(xs * arc_weight, ys * arc_weight, np.transpose(intensity_cart), cmap = cmap)
    result.set_clim(clim[0], clim[1])

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad) * arc_weight, color = "black")
    ax.add_artist(circle)

    # Add beam size
    beam = plot.Circle((-2 * arc_weight, -2 * arc_weight), (beam_size / 2) * arc_weight, color = "white")
    fig.gca().add_artist(beam)

    # Add planet orbit
    planet_orbit = plot.Circle((0, 0), 1 * arc_weight, color = "white", fill = False, alpha = 0.8, linestyle = "dashed", zorder = 50)
    ax.add_artist(planet_orbit)

    # Locate Planet
    if shift < -len(dust_theta):
        shift += len(dust_theta)
    planet_theta = dust_theta[shift]
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
    plot.scatter(planet_x * arc_weight, planet_y * arc_weight, c = "white", s = int(70 * planet_size), marker = "D", zorder = 100) # planet

    # Annotate Axes
    if arc:
        unit = "^{\prime\prime}"
    else:
        unit = "r_\mathrm{p}"

    ax.set_xlabel(r"$x$ [$%s$]" % unit, fontsize = fontsize)
    if frame_i == 1:
        ax.set_ylabel(r"$y$ [$%s$]" % unit, fontsize = fontsize)

    # Axes
    box_size = args.box * arc_weight
    ax.set_xlim(-box_size, box_size)
    ax.set_ylim(-box_size, box_size)
    ax.set_aspect('equal')

    ax.spines['bottom'].set_color('w'); ax.spines['top'].set_color('w'); ax.spines['left'].set_color('w'); ax.spines['right'].set_color('w')
    ax.tick_params(colors = 'white', labelcolor = 'black', width = 1, length = 5)

    # Label
    taper_title = r"$T_\mathrm{growth} = %d$" % taper_time
    plot.text(0.0 * box_size, 2 * arc_weight, taper_title, fontsize = fontsize, color = 'white', horizontalalignment = 'center', bbox=dict(facecolor = 'black', edgecolor = 'white', pad = 10.0))

    # Title
    title = r"$t$ $=$ $%.1f$  [$m_p(t)$ $=$ $%.2f$ $M_J$]" % (orbit, current_mass)
    plot.title("%s" % (title), y = 1.015, fontsize = fontsize + 1)

    # Super Title
    left_x = -1.2 * box_size; line_y = 1.24 * box_size; linebreak = 0.2 * box_size
    right_x = 1.2 * box_size
    if frame_i == 1:
        #line1 = r'$M_p = %d$ $M_J$' % planet_mass
        #line2 = r'$\nu = 10^{%d}$' % round(np.log(viscosity) / np.log(10), 0)
        #plot.text(left_x, line_y + 1.2 * linebreak, line1, horizontalalignment = 'left', fontsize = fontsize + 2)
        #plot.text(left_x, line_y + 0.2 * linebreak, line2, horizontalalignment = 'left', fontsize = fontsize + 2)
        line1 = r'$M_p = %d$ $M_J$' % planet_mass
        plot.text(left_x, line_y + 0.2 * linebreak, line1, horizontalalignment = 'left', fontsize = fontsize + 2)
    elif frame_i == 2:
        #line3 = r"$%.02f^{\prime\prime} \times \ \ %.02f^{\prime\prime}$" % (arc_beam, arc_beam)
        #plot.text(right_x, line_y + 0.7 * linebreak, line3, horizontalalignment = 'right', fontsize = fontsize + 2)
        line2 = r'$\nu = 10^{%d}$' % round(np.log(viscosity) / np.log(10), 0)
        plot.text(right_x, line_y + 0.2 * linebreak, line2, horizontalalignment = 'right', fontsize = fontsize + 2)

    # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
    if colorbar:
        # Only for last frame
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "8%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)
        if frame_i % 2 == 0:
            cbar.set_label("Normalized Intensity", fontsize = fontsize, rotation = 270, labelpad = 25)

        #if frame_i != num_frames:
        #    fig.delaxes(cax) # to balance out frames that don't have colorbar with the one that does

    # Return to previous directory
    os.chdir(cwd)

    return ax
    
def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (12, 6), dpi = dpi)
    gs = gridspec.GridSpec(1, 2)

    frame_str = ""
    for i, frame_i in enumerate(frame_range):
        ax = fig.add_subplot(gs[i])
        ax = add_to_plot(frame_i, fig, ax, len(frame_range), i + 1)
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    #### Finish Plot ####
    title = r"$%.02f^{\prime\prime} \times \ \ %.02f^{\prime\prime}$" % (arc_beam, arc_beam)
    #title = r"$\mathrm{Beam:\ }\ %.02f^{\prime\prime} \times \ \ %.02f^{\prime\prime}$" % (arc_beam, arc_beam)
    #title = r"$N_\mathrm{r} \times \ N_\mathrm{\phi} = %d \times \ %d$" % (old_num_rad, old_num_theta)
    fig.suptitle(title, y = 0.97, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)

    # Save and Close
    plot.tight_layout()

    # Save, Show,  and Close
    if version is None:
        save_fn = "%s/intensityCartGridTaperComparison_%s.png" % (save_directory, frame_str)
    else:
        save_fn = "%s/v%04d_intensityCartGridTaperComparison_%s.png" % (save_directory, version, frame_str)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)