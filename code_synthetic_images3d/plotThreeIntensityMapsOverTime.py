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
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pylab import rcParams
from pylab import fromfile

import util
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot convolved intensity maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = 3,
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "threeIntensityMapsInTime",
                         help = 'save directory (default: threeIntensityMapsInTime)')

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

    parser.add_argument('--cbar', dest = "colorbar", action = 'store_true', default = False,
                         help = 'include colorbar (default: no colorbar)')
    parser.add_argument('--sup', dest = "supertitle", action = 'store_true', default = False,
                         help = 'include super title (default: do not)')
    parser.add_argument('--title', dest = "optional_title", default = None,
                         help = 'optional title (default: None)')

    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "inferno",
                         help = 'color map (default: inferno)')
    parser.add_argument('--cmax', dest = "cmax", type = int, default = None,
                         help = 'maximum density in colorbar (default: 2.5)')

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
fn = "id%04d_par.p" % args.id_number
fargo_par = pickle.load(open(fn, "rb"))

p = fargo_par["p"]

num_rad = p.ny; num_theta = p.nx
r_min = p.ymin; r_max = p.ymax

surface_density_zero = p.sigma0
dust_surface_density_zero = p.sigma0 * p.epsilon

scale_height = p.aspectratio
viscosity = p.nu

dt = p.ninterm * p.dt

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

planet_radius = fargo_par["Radius"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

real_fargo_par = util.get_pickled_parameters(directory = "../../")
jupiter_mass = 1e-3
planet_mass = real_fargo_par["PlanetMass"] / jupiter_mass
accretion = real_fargo_par["Accretion"]
taper_time = p.masstaper

### Get Input Parameters ###

# Frames
frame_range = args.frames

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

# Planet File
# Data
data = np.loadtxt("../../planet0.dat")
times = data[:, 0]; base_mass = data[:, 7]
accreted_mass = data[:, 8] / jupiter_mass

##### PLOTTING #####

def make_plot(frames, show = False):
    # Set up figure
    fig = plot.figure(figsize = (12, 4), dpi = dpi)

    # Indvidual Subplots
    def add_to_plot(i):
        # Identify Subplot
        frame = frames[i]; number = i + 1
        ax = plot.subplot(1, 3, number)

        # Data
        intensity_cart = util.read_data(frame, 'cartesian_intensity', fargo_par, id_number = id_number)
        xs, ys, xs_grid, ys_grid = sq.get_cartesian_grid(rad)

        # Get Shift
        gas_density = util.read_data(frame, 'gas', fargo_par, directory = "../../")

        # Shift gas density with center of dust density
        shift = az.shift_away_from_minimum(gas_density, fargo_par)

        # Normalize
        if normalize:
            intensity_cart /= np.max(intensity_cart)

        # Arcseconds or Planet Radii
        if arc:
            arc_weight = planet_radius / distance # related to parallax
        else:
            arc_weight = 1

        ### Plot ###
        result = plot.pcolormesh(xs * arc_weight, ys * arc_weight, np.transpose(intensity_cart), cmap = cmap)
        result.set_clim(clim[0], clim[1])

        # Get rid of interior
        circle = plot.Circle((0, 0), min(rad) * arc_weight, color = "black")
        ax.add_artist(circle)

        # Add beam size
        beam = plot.Circle((1.7 * arc_weight, 1.7 * arc_weight), (beam_size / 2) * arc_weight, color = "white")
        fig.gca().add_artist(beam)

        # Add planet orbit
        planet_orbit = plot.Circle((0, 0), 1 * arc_weight, color = "white", fill = False, alpha = 0.8, linestyle = "dashed", zorder = 50)
        ax.add_artist(planet_orbit)

        # Locate Planet
        if shift < -len(theta):
            shift += len(theta)
        planet_theta = theta[shift]
        planet_theta += (np.pi / 2.0) # Note: the conversion from polar to cartesian rotates everything forward by 90 degrees
        planet_theta = planet_theta % (2 * np.pi) - (np.pi) # Keep -np.pi < theta < np.pi

        planet_x = np.cos(planet_theta)
        planet_y = np.sin(planet_theta)

        # Label star and planet
        time = fargo_par["Ninterm"] * fargo_par["DT"]
        orbit = (time / (2 * np.pi)) * frame
        if orbit >= taper_time:
            current_mass = planet_mass
        else:
            current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

        current_mass += accreted_mass[frame]

        plot.scatter(0, 0, c = "white", s = 250, marker = "*", zorder = 100) # star
        plot.scatter(planet_x * arc_weight, planet_y * arc_weight, c = "white", s = 60, marker = "D", zorder = 100) # planet

        # Axes
        box_size = args.box * arc_weight

        ax.set_xlim(-box_size, box_size)
        ax.set_ylim(-box_size, box_size)
        ax.set_aspect('equal')

        ax.spines['bottom'].set_color('w'); ax.spines['top'].set_color('w'); ax.spines['left'].set_color('w'); ax.spines['right'].set_color('w')
        ax.tick_params(colors = 'white', labelcolor = 'black', width = 1, length = 5, direction = "in")

        # Annotate Axes
        if arc:
            ticks = np.arange(-0.3, 0.31, 0.1)
            ax.set_xticks(ticks); ax.set_xticklabels(["", "-0.2", "", "0.0", "", "0.2", ""])
            ax.set_yticks(ticks)
            unit = "^{\prime\prime}"
        else:
            unit = "r_\mathrm{p}"

        ax.set_xlabel(r"$x$ [$%s$]" % unit, fontsize = fontsize)
        if number == 1:
            ax.set_ylabel(r"$y$ [$%s$]" % unit, fontsize = fontsize)

        # Title
        x_min = plot.xlim()[0]; x_max = plot.xlim()[-1]
        x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
        x_shift = 0.35; extra = 0.17
        y_text = 1.16; y_shift = 0.10

        alpha_coefficent = "3"
        if scale_height == 0.08:
            alpha_coefficent = "1.5"
        elif scale_height == 0.04:
            alpha_coefficent = "6"

        if i == 0:
            text1 = r"$h = %.2f$" % (scale_height)
            plot.text(x_min - x_shift * x_range, (y_text + y_shift) * plot.ylim()[-1], text1, horizontalalignment = 'left', fontsize = fontsize + 1)
            text2 = r"$\alpha \approx %s \times 10^{%d}$" % (alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2)
            plot.text(x_min - x_shift * x_range, (y_text) * plot.ylim()[-1], text2, horizontalalignment = 'left', fontsize = fontsize + 1)
        if i == 2:
            text3 = args.optional_title
            #text4 = r"%d$\times$%d (2-D)" % (num_rad, num_theta)

            if text3 is not None:
              #plot.text(x_max + x_shift * x_range, (y_text + 0.5 * y_shift) * plot.ylim()[-1], text3, horizontalalignment = 'right', fontsize = fontsize + 1)
              plot.text(x_max + (x_shift + extra) * x_range, (y_text + y_shift + 0.01) * plot.ylim()[-1], text3, horizontalalignment = 'right', fontsize = fontsize + 1)
              #plot.text(x_max + (x_shift + extra) * x_range, (y_text + 0.01) * plot.ylim()[-1], text4, horizontalalignment = 'right', fontsize = fontsize + 1)

        title = r"$t = %d$ [$m_\mathrm{p}=%.2f$ $M_\mathrm{J}$]" % (orbit, current_mass)
        plot.title("%s" % (title), y = 1.035, fontsize = fontsize)

        # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
        if args.colorbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size = "6%", pad = 0.2)
            #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
            cbar = fig.colorbar(result, cax = cax)
            cbar.set_label(r"Normalized Intensity", fontsize = fontsize, rotation = 270, labelpad = 25)

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
    #title = r"$h = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)

    #beam_diameter = fargo_par["Beam"] * fargo_par["Radius"] / fargo_par["Distance"]
    if args.supertitle:
        #title = r'$h = %.2f$   $\Sigma = %.3e$  (2-D)  [$%.3f^{\prime\prime}$]' % (scale_height, fargo_par["p"].sigma0, arc_beam)
        title = r"$\Sigma_0$ $/$ $\Sigma_\mathrm{base} = %.1f$    $M_\mathrm{p} = %.2f$ $M_\mathrm{Jup}$    $%.3f^{\prime\prime}$" % (surface_density_zero / surface_density_base, final_planet_mass, arc_beam)
        #plot.suptitle("%s" % (title), y = 1.15, fontsize = fontsize + 2, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))
        plot.suptitle("%s" % (title), y = 1.32, fontsize = fontsize + 2, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))

    # Tighten!
    plot.tight_layout()

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/id%04d_intensityCartGrid_%04d-%04d-%04d.png" % (save_directory, id_number, frames[0], frames[1], frames[2])
    else:
        save_fn = "%s/v%04d_id%04d_intensityCartGrid_%04d-%04d-%04d.png" % (save_directory, version, id_number, frames[0], frames[1], frames[2])
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi, pad_inches = 0.15)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####
make_plot(frame_range, show = show)