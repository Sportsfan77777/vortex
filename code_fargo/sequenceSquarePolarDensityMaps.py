"""
plot 2-D density maps

python plotDensityMaps.py
python plotDensityMaps.py frame_number
python plotDensityMaps.py -1 <<<===== Plots a sample
python plotDensityMaps.py -m
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np

from scipy.interpolate import interp1d as interpolate
from scipy.ndimage import map_coordinates

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pylab import rcParams
from pylab import fromfile

import util
from readTitle import readTitle

save_directory = "squareDensityMapSequences"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

planet_mass = float(fargo_par["PlanetMass"])
taper_time = float(fargo_par["MassTaper"])

### Converter ###

def polar_to_cartesian(data, rs, thetas, order = 3):
    # Source: http://stackoverflow.com/questions/2164570/reprojecting-polar-to-cartesian-grid
    # Set up xy-grid
    max_r = rs[-1]
    resolution = 2 * len(rad)

    xs = np.linspace(-max_r, max_r, resolution)
    ys = np.linspace(-max_r, max_r, resolution)

    xs_grid, ys_grid = np.meshgrid(xs, ys)

    # Interpolate rt-grid

    interpolated_rs = interpolate(rs, np.arange(len(rs)), bounds_error = False)
    interpolated_thetas = interpolate(thetas, np.arange(len(thetas)), bounds_error = False)

    # Match up xy-grid with rt-grid

    new_rs = np.sqrt(np.power(xs_grid, 2) + np.power(ys_grid, 2))
    new_thetas = np.arctan2(ys_grid, xs_grid)

    # Convert from [-pi, pi] to [0, 2 * pi]
    new_thetas[new_thetas < 0] = new_thetas[new_thetas < 0] + 2 * np.pi

    new_interpolated_rs = interpolated_rs(new_rs.ravel())
    new_interpolated_thetas = interpolated_thetas(new_thetas.ravel())

    # Fix Bounds (outside of polar image, but inside cartesian image)
    new_interpolated_rs[new_rs.ravel() > max(rs)] = len(rs) - 1
    new_interpolated_rs[new_rs.ravel() < min(rs)] = 0

    cart_data = map_coordinates(data, np.array([new_interpolated_rs, new_interpolated_thetas]), order = order).reshape(new_rs.shape)

    return xs_grid, ys_grid, cart_data

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "RdYlBu_r"
clim = [0, 2]

fontsize = 14
my_dpi = 100

def add_to_plot(ax, frame, num_frames, frame_i):
    print frame, num_frames, frame_i

    # Declare Subplot
    #ax = plot.subplot(1, num_frames, frame_i, sharex = prev_ax, sharey = prev_ax, aspect = "equal")

    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Mass
    if orbit >= taper_time:
        current_mass = planet_mass / 0.001
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (orbit / taper_time)), 2) * (planet_mass / 0.001)

    # Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero
    xs_grid, ys_grid, density_cart = polar_to_cartesian(density, rad, theta)

    # Axes
    sq = 2.5
    ax.set_xlim(-sq, sq)
    ax.set_ylim(-sq, sq)
    #plot.axes().set_aspect('equal')

    ### Plot ###
    result = plot.pcolormesh(xs_grid, ys_grid, np.transpose(density_cart), cmap = cmap)
    result.set_clim(clim[0], clim[1])

    # Get rid of interior
    circle = plot.Circle((0, 0), min(rad), color = "black")
    ax.add_artist(circle)

    # Add minor grid lines
    alpha = 0.25
    dashes = [1, 5]
    ax.grid(b = True, which = 'major', color = "black", dashes = dashes, alpha = alpha)
    ax.grid(b = True, which = 'minor', color = "black", dashes = dashes, alpha = alpha)
    ax.minorticks_on()

    # Annotate
    #this_title = readTitle()
    #plot.xlabel("x", fontsize = fontsize)
    #plot.ylabel("y", fontsize = fontsize)
    ax.set_title(r"$t = %d$, $m_p(t) = %.2f$ $M_J$" % (orbit, current_mass), fontsize = fontsize + 1)

    # Add Colorbar
    if frame_i == num_frames:
        # Only for last frame
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])

        fig.colorbar(result)

    return ax
    

def finish_plot(frame_range, show = True):
    frame_str = ""
    for i, frame in enumerate(frame_range):
        if i == 0:
            frame_str += "%d" % frame
        else:
            frame_str += "-%d" % frame

    # Save and Close
    plot.savefig("%s/densityMapSequence_%s.png" % (save_directory, frame_str), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


##### Plot Files #####

widths = [700, 1100, 1500, 1900]

if len(sys.argv) > 1:
    frame_range = [int(frame) for frame in sys.argv[1:]]

    # Set up figure
    #fig = plot.figure(figsize = (widths[len(frame_range) - 1] / my_dpi, 600 / my_dpi), dpi = my_dpi)
    fig = plot.figure(dpi = my_dpi)
    gs = gridspec.GridSpec(1, len(frame_range))

    for i, frame in enumerate(frame_range):
        if i == 0:
            ax = fig.add_subplot(gs[0])
        else:
            ax = fig.add_subplot(gs[i], sharey = True)

        add_to_plot(ax, frame, len(frame_range), i + 1)

    finish_plot(frame_range)

