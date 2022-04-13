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

from pylab import rcParams
from pylab import fromfile

import util
import utilVorticity
import azimuthal as az
#from readTitle import readTitle

from advanced import Parameters
from reader_mpiio import Fields

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
    parser.add_argument('--dir', dest = "save_directory", default = "meridionalVorticityMaps",
                         help = 'save directory (default: meridionalVorticityMaps)')
    parser.add_argument('-m', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpio output files (default: use dat)')
    parser.add_argument('--merge', dest = "merge", type = int, default = 0,
                         help = 'number of cores needed to merge data outputs (default: 0)')

    # Quantity to plot
    parser.add_argument('--rossby', dest = "rossby", action = 'store_true', default = False,
                         help = 'plot rossby number instead of vorticity (default: plot vorticity)')
    parser.add_argument('--residual', dest = "residual", action = 'store_false', default = True,
                         help = 'use v_theta or v_theta - v_kep (default: use residual)')

    parser.add_argument('-z', dest = "sliver", type = int, default = 0,
                         help = 'sliver above midplane to plot (default: 0)')

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
    parser.add_argument('--scale', dest = "quiver_scale", type = float, default = 0.25,
                         help = 'bigger scale means smaller arrow (default: 1)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--cmap', dest = "cmap", default = "magma",
                         help = 'color map (default: magma)')
    parser.add_argument('--crange', dest = "c_lim", type = float, nargs = 2, default = None,
                         help = 'range in colorbar (default: [-0.2, 0])')

    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
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
viscosity = 1e-7 #p.nu

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
frame_range = util.get_frame_range(args.frames)

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
if args.c_lim is None:
    clim = [-0.2, 0]
else:
    clim = [args.c_lim[0], args.c_lim[1]]

fontsize = args.fontsize
dpi = args.dpi

# Planet File
# Data
data = np.loadtxt("planet0.dat")
times = data[:, 0]; base_mass = data[:, 7]
accreted_mass = data[:, 8] / jupiter_mass

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
    shifted_vorticity = np.roll(vorticity, shift_c, axis = -1)
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

def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    if mpi:
      density = Fields("./", 'gas', frame).get_field("dens").reshape(num_z, num_rad, num_theta)
      vrad = Fields("./", 'gas', frame).get_field("vy").reshape(num_z, num_rad, num_theta)
      vtheta = Fields("./", 'gas', frame).get_field("vx").reshape(num_z, num_rad, num_theta)
    else:
      density = fromfile("gasdens%d.dat" % frame).reshape(num_z, num_rad, num_theta)
      vrad = (fromfile("gasvy%d.dat" % frame).reshape(num_z, num_rad, num_theta)) # add a read_vrad to util.py!
      vtheta = (fromfile("gasvx%d.dat" % frame).reshape(num_z, num_rad, num_theta)) # add a read_vrad to util.py!
    midplane_density = density[num_z / 2 + args.sliver, :, :]
    midplane_vrad = vrad[num_z / 2 + args.sliver, :, :]
    midplane_vtheta = vtheta[num_z / 2 + args.sliver, :, :]

    dz = z_angles[1] - z_angles[0]
    surface_density = np.sum(density[:, :, :], axis = 0) * dz

    normalized_midplane_density = midplane_density / surface_density_zero # / np.sqrt(2.0 * np.pi) / scale_height_function[:, None]
    normalized_density = surface_density / surface_density_zero # / np.sqrt(2.0 * np.pi) / scale_height_function[:, None]

    #vorticity = utilVorticity.velocity_curl(midplane_vrad, midplane_vtheta, rad, theta, rossby = rossby, residual = residual)
    vorticity = utilVorticity.velocity_curl(vrad, vtheta, rad, theta, rossby = rossby, residual = residual)
    voriticity_cut = vorticity[:, :, num_theta / 2.0]

    # Shift
    if center:
       normalized_density, vorticity, shift_c = shift_density(normalized_density, vorticity, fargo_par, reference_density = normalized_density)

    ### Plot ###
    x = rad
    y = z_angles # * (180.0 / np.pi)
    result = ax.pcolormesh(x, y, np.transpose(voriticity_cut), cmap = cmap)

    cbar = fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    if use_contours:
        levels = np.linspace(low_contour, high_contour, num_levels)
        colors = generate_colors(num_levels)
        plot.contour(x, y, np.transpose(normalized_density), levels = levels, origin = 'upper', linewidths = 1, colors = colors, alpha = 0.8)

    if quiver:
        # Velocity
        radial_velocity = vrad
        azimuthal_velocity = vtheta
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

        plot.quiver(x_q[::rate_x], y_q[::rate_y], u[::rate_y,::rate_x], v[::rate_y,::rate_x], scale = scale)

    # Axes
    plot.xlim(x_min, x_max)
    plot.ylim(z_min, z_max)

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    current_mass += accreted_mass[frame]

    #title = readTitle()

    unit = "r_\mathrm{p}"
    plot.xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
    plot.ylabel(r"$\theta$", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)
    #plot.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    #text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    #text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')

    # Label colorbar
    if rossby:
       cbar_name = r"$\mathrm{Rossby}$ $\mathrm{number}$"
    else:
       cbar_name = r"$\mathrm{Vorticity}$"
    cbar.set_label(cbar_name, fontsize = fontsize, rotation = 270, labelpad = 25)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/vorticityMap_%04d.png" % (save_directory, frame)
    else:
        save_fn = "%s/v%04d_vorticityMap_%04d.png" % (save_directory, version, frame)
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

