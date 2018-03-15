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
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "cartesianIntensityMaps",
                         help = 'save directory (default: cartesianIntensityMaps)')

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

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

planet_radius = fargo_par["Radius"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

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

###############################################################################

##### PLOTTING #####

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    intensity_cart = util.read_data(frame, 'cartesian_intensity', fargo_par, id_number = id_number)
    xs, ys, xs_grid, ys_grid = sq.get_cartesian_grid(rad)

    # Get Shift
    dust_fargo_par = util.get_pickled_parameters(directory = "../../../cm-size") ## shorten name?
    ######## Need to extract parameters, and add 'rad' and 'theta' ########
    dust_rad = np.linspace(dust_fargo_par['Rmin'], dust_fargo_par['Rmax'], dust_fargo_par['Nrad'])
    dust_theta = np.linspace(0, 2 * np.pi, dust_fargo_par['Nsec'])
    dust_fargo_par['rad'] = dust_rad; dust_fargo_par['theta'] = dust_theta
    gas_surface_density_zero = dust_fargo_par['Sigma0']

    dust_density = util.read_data(frame, 'dust', dust_fargo_par, id_number = id_number, directory = "../../../cm-size")

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

    # Axes
    box_size = args.box * arc_weight
    plot.xlim(-box_size, box_size)
    plot.ylim(-box_size, box_size)
    plot.axes().set_aspect('equal')

    ax.spines['bottom'].set_color('w'); ax.spines['top'].set_color('w'); ax.spines['left'].set_color('w'); ax.spines['right'].set_color('w')
    ax.tick_params(colors = 'white', labelcolor = 'black', width = 1, length = 5)

    # Annotate Axes
    if arc:
        unit = "^{\prime\prime}"
    else:
        unit = "r_\mathrm{p}"

    ax.set_xlabel(r"$x$ [$%s$]" % unit, fontsize = fontsize)
    ax.set_ylabel(r"$y$ [$%s$]" % unit, fontsize = fontsize)

    # Title
    title = r"$t$ $=$ $%.1f$  [$m_p(t)$ $=$ $%.2f$ $M_J$]" % (orbit, current_mass)
    plot.title("%s" % (title), y = 1.015, fontsize = fontsize + 1)

    #plot.title("Intensity Map (t = %.1f)" % (orbit), fontsize = fontsize + 1)

    taper_title = r"$T_\mathrm{growth} = %d$" % taper_time
    plot.text(0.0 * box_size, 2 * arc_weight, taper_title, fontsize = fontsize, color = 'white', horizontalalignment = 'center', bbox=dict(facecolor = 'black', edgecolor = 'white', pad = 10.0), zorder = 100)

    # Add Colorbar (Source: http://stackoverflow.com/questions/23270445/adding-a-colorbar-to-two-subplots-with-equal-aspect-ratios)
    if colorbar:
        # Only for last frame
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "8%", pad = 0.2)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(result, cax = cax)
        cbar.set_label("Normalized Intensity", fontsize = fontsize, rotation = 270, labelpad = 25)

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

def old_make_plot(frame, show = False):
    """
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * i

    # Mass
    if orbit >= taper_time:
        current_mass = planet_mass / 0.001
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * (planet_mass / 0.001)

    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    ax = fig.add_subplot(111)

    # Data  <<<---- Read in data outside of this function.
    data = np.loadtxt("intensitymap.out")
    intensity = data[:, -1].reshape(num_rad, num_theta)
    intensity = clear_inner_disk(intensity)
    xs, ys, xs_grid, ys_grid, intensity_cart = polar_to_cartesian(intensity, rad, theta)

    # Convolve Data
    convolved_intensity = convolve_intensity(intensity_cart)
    contrast, maximum, opposite = record_contrast(convolved_intensity, xs, ys)

    if axis == "zoom":
        prefix = "zoom_"
        sq = 2.5
    else:
        prefix = ""
        sq = np.max(xs_grid)

    plot.xlim(-sq, sq)
    plot.ylim(-sq, sq)
    plot.axes().set_aspect('equal')

    ### Plot ###
    result = ax.pcolormesh(xs_grid, ys_grid, np.transpose(convolved_intensity), cmap = cmap)
    cbar = fig.colorbar(result)

    clim = [np.percentile(intensity[clim_in : clim_out, :], clim_low), np.percentile(intensity[clim_in : clim_out, :], clim_high)]
    result.set_clim(clim[0], clim[1])

    cbar.set_label(r"Intensity", fontsize = fontsize + 2, rotation = 270, labelpad = 20)

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    fig.gca().add_artist(circle)

    # Add beam size
    beam = plot.Circle((-2, -2), beam_size, color = "white")
    fig.gca().add_artist(beam)

    # Label star and planet
    planet_size = (current_mass / (planet_mass / 0.001))
    plot.scatter(0, 0, c = "white", s = 300, marker = "*", zorder = 100) # star
    plot.scatter(0, 1, c = "white", s = int(70 * planet_size), marker = "D") # planet

    # Add minor grid lines
    alpha = 0.2
    dashes = [1, 5]
    #plot.grid(b = True, which = 'major', color = "black", dashes = dashes, alpha = alpha)
    #plot.grid(b = True, which = 'minor', color = "black", dashes = dashes, alpha = alpha)
    plot.minorticks_on()

    # Annotate
    title1 = r"$m_p = %d$ $M_{Jup}$, $\nu_{disk} = 10^{%d}$, $T_{growth} = %d$ $\rm{orbits}$" % (int(planet_mass / 0.001), int(np.log(viscosity) / np.log(10)), taper_time)
    title2 = r"$t = %d$ $\rm{orbits}}$, $m_p(t) = %.2f$ $M_{Jup}$ | $r = %d$ $\rm{AU}$, $\lambda = %d$ $\rm{\mu m}$" % (orbit, current_mass, int(radius), wavelength)
    #plot.xlabel("x", fontsize = fontsize)
    #plot.ylabel("y", fontsize = fontsize)
    plot.title("%s" % (title2), y = 1.01, fontsize = fontsize)
    plot.text(0.0, 3.14, title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Write contrast
    plot.text(-2, 2, "Contrast: %.2f = %.2e / %.2e" % (contrast, maximum, opposite), color = "white", fontsize = fontsize)

    # Save and Close
    plot.savefig("%sconvolvedIntensityMap%04d_r%d_at%d.png" % (prefix, i, int(radius), wavelength), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)
    """


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



