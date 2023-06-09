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
import scipy

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
    parser.add_argument('--dir', dest = "save_directory", default = "timeAveragedVerticalAlpha",
                         help = 'save directory (default: timeAveragedVerticalAlpha)')
    parser.add_argument('-m', dest = "mpi", action = 'store_true', default = False,
                         help = 'use .mpiio output files (default: do not use mpi)')
    parser.add_argument('--merge', dest = "merge", type = int, default = 0,
                         help = 'number of cores needed to merge data outputs (default: 0)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')
    parser.add_argument('-w', dest = "window", type = int, default = 31,
                         help = 'window size for smoothing (default: 31)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--y_range', dest = "y_range", type = float, nargs = 2, default = [-0.25, 0],
                         help = 'range in y-axis (default: [-0.5, 0.0])')

    parser.add_argument('--compare', dest = "compare", default = None,
                         help = 'compare to another directory (default: do not do it!)')

    # Quantity to plot for maximum condition
    parser.add_argument('--rossby', dest = "rossby", action = 'store_false', default = True,
                         help = 'plot rossby number instead of vorticity (default: plot vorticity)')
    parser.add_argument('--residual', dest = "residual", action = 'store_false', default = True,
                         help = 'use v_theta or v_theta - v_kep (default: do not use residual)')

    parser.add_argument('-z', dest = "sliver", type = int, default = 0,
                         help = 'sliver above midplane to plot (default: 0)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 21,
                         help = 'fontsize of plot annotations (default: 21)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 4,
                         help = 'fontsize of plot annotations (default: 3)')
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

planet_mass = 1.0
taper_time = p.masstaper

scale_height = p.aspectratio
flaring_index = p.flaringindex
viscosity = 1e-7 #p.nu

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

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
z_angles = np.linspace(z_min, z_max, num_z)

version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
#max_y = args.max_y
y_range = args.y_range

# Quantity to Plot
rossby = args.rossby
residual = args.residual

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

# Planet File
# Data
data = np.loadtxt("planet0.dat")
times = data[:, 0]; base_mass = data[:, 7]
accreted_mass = data[:, 8] / jupiter_mass

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

def generate_colors(n):
    c = ['yellow', 'b', 'firebrick', 'w', 'green']
    colors = []
    for i in range(n):
        colors.append(c[i % len(c)])
    return colors

##### PLOTTING #####

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    total_averaged_reynolds_stress_z_theta = np.zeros((num_z, num_rad))
    rolling_averaged_reynolds_stress_z_theta = np.zeros((num_z, num_rad))

    # Time-averaged Data

    time_directory = "timeAverages"
    time_density = np.load("%s/time_averaged_density_%04d-%04d-%04d.npy" % (time_directory, args.frames[0], args.frames[1], args.frames[2]))
    time_vz = np.load("%s/time_averaged_vz_%04d-%04d-%04d.npy" % (time_directory, args.frames[0], args.frames[1], args.frames[2]))
    time_vrad = np.load("%s/time_averaged_vy_%04d-%04d-%04d.npy" % (time_directory, args.frames[0], args.frames[1], args.frames[2]))
    time_vtheta = np.load("%s/time_averaged_vx_%04d-%04d-%04d.npy" % (time_directory, args.frames[0], args.frames[1], args.frames[2])) 

    for frame in frame_range:
        print frame

        # Data
        density = fromfile("gasdens%d.dat" % frame).reshape(num_z, num_rad, num_theta)
        vz = (fromfile("gasvz%d.dat" % frame).reshape(num_z, num_rad, num_theta)) # add a read_vrad to util.py!
        vrad = (fromfile("gasvy%d.dat" % frame).reshape(num_z, num_rad, num_theta)) # add a read_vrad to util.py!
        vtheta = (fromfile("gasvx%d.dat" % frame).reshape(num_z, num_rad, num_theta)) # add a read_vrad to util.py!

        # Reynolds Stress
        vz_component = vz - (np.average(time_vz, axis = -1))[:, :, None]
        vtheta_component = vtheta - (np.average(time_vtheta, axis = -1))[:, :, None]

        reynolds_stress_z_theta = density * vz_component * vtheta_component
        averaged_reynolds_stress_z_theta = np.average(reynolds_stress_z_theta, axis = -1)

        # Sum over time
        total_averaged_reynolds_stress_z_theta += averaged_reynolds_stress_z_theta

    # Time Average
    total_averaged_reynolds_stress_z_theta /= len(frame_range)

    # Rolling Average
    dr = rad[1] - rad[0]
    window_size = args.window # scale_height / dr
    poly_order = 1

    smoothed_reynolds_stress = scipy.signal.savgol_filter(total_averaged_reynolds_stress_z_theta, window_size, poly_order)

    # Add up columns
    dz = z_angles[1] - z_angles[0]
    #surface_smoothed_reynolds_stress = np.sum(smoothed_reynolds_stress[:, :], axis = 0) * dz

    # Calculate alpha
    q = 1.0 # Power of temperature power law

    omega = np.power(rad, -1.5)
    omega_gradient = -1.5 * np.power(rad, -2.5)
    scale_height_function = scale_height * rad
    sound_speed_function = scale_height_function * omega
    zs = z_angles - (np.pi / 2.0)

    surface_density = np.average(np.sum(density[:, :, :], axis = 0) * dz, axis = -1)

    coefficient = 0.5 * zs[:, None] * q * scale_height / scale_height_function
    averaged_pressure = np.average(time_density, axis = -1) * np.power(sound_speed_function, 2.0)

    measured_alpha = smoothed_reynolds_stress / averaged_pressure / coefficient

    r_choice = 1.4
    r_i = np.searchsorted(rad, r_choice)

    ### Plot ###
    x = zs / scale_height
    y = measured_alpha[:, r_i]

    result, = plot.plot(x, y, linewidth = linewidth, c = "b", label = "min", zorder = 90)
    result2, = plot.plot(x, -y, linewidth = linewidth, c = "r", label = "min", zorder = 10, alpha = 0.7)

    # Axes
    plot.xlim(x[0], x[-1])
    plot.ylim(10**(-2), 10**(0))
    plot.yscale('log')
    #plot.yticks(np.arange(y_range[0], y_range[1] + 1e-9, 0.005))

    # Annotate Axes
    orbit = (dt / (2 * np.pi)) * frame

    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass

    current_mass += accreted_mass[frame]

    #title = readTitle()

    unit = "r_\mathrm{p}"
    ax.set_xlabel(r"z / H [$%s$]" % unit, fontsize = fontsize)
    ax.set_ylabel(r"$\alpha_z \approx <\Delta v_z \Delta v_{\phi}> / c_s^2$", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)
    title1 = r"$h/r = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    title2 = r"$t = %d$ $\mathrm{orbits}}$  [$m_\mathrm{p}(t)\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)
    #ax.text(x_mid, y_text * plot.ylim()[0], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')

    # Save, Show, and Close
    directory_name = os.getcwd().split("/")[-1]

    if version is None:
        save_fn = "%s/%s_timeAveragedVerticalAlpha_%04d-%04d-%04d.png" % (save_directory, directory_name, args.frames[0], args.frames[1], args.frames[2])
    else:
        save_fn = "%s/v%04d_%s_timeAveragedVerticalAlpha_%04d-%04d-%04d.png" % (save_directory, version, directory_name, args.frames[0], args.frames[1], args.frames[2])
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plots! #####

make_plot(show = show)
