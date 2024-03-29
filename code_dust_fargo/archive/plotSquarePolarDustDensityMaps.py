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

from pylab import rcParams
from pylab import fromfile

import util
from readTitle import readTitle

import colormaps as cmaps
plot.register_cmap(name = 'viridis', cmap = cmaps.viridis)
plot.register_cmap(name = 'inferno', cmap = cmaps.inferno)
plot.register_cmap(name = 'plasma', cmap = cmaps.plasma)
plot.register_cmap(name = 'magma', cmap = cmaps.magma)

save_directory = "squarePolarDustDensityMaps"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = float(fargo_par["Nrad"])
num_theta = float(fargo_par["Nsec"])

rad = np.linspace(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]), num_rad + 1)
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"]) / 100.0
scale_height = float(fargo_par["AspectRatio"])
viscosity = float(fargo_par["Viscosity"])

planet_mass = float(fargo_par["PlanetMass"])
taper_time = int(float(fargo_par["MassTaper"]))

# Size Labels
size = sys.argv[2]
size_label = {}
size_label["cm"] = r"$\rm{cm}$"
size_label["hcm"] = r"$3$ " + r"$\rm{mm}$"
size_label["mm"] = r"$\rm{mm}$"
size_label["hmm"] = r"$0.3$ "r"$\rm{mm}$"
size_label["hum"] = r"$100$ " + r"$\rm{\mu m}$"
size_label["um"] = r"$\rm{\mu m}$"

### Converter ###

def clear_inner_disk(data):
    # get rid of inner disk (r < outer_limit)
    filtered_data = np.copy(data)

    outer_limit = np.searchsorted(rad, 1.05)
    filtered_data[:outer_limit] = 0

    return filtered_data

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
cmap = "inferno"
clim = [0, 5]

fontsize = 15
my_dpi = 100

def make_plot(frame, size, show = False):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    print frame
    def choose_axis(i, axis):
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

        # Data
        if len(sys.argv) > 3:
            density = pickle.load(open("shifted_gasddens%d_%s.p" % (i, size), 'r')) / surface_density_zero
        else:
            density = (fromfile("gasddens%d_%s.dat" % (i, size)).reshape(num_rad, num_theta)) / surface_density_zero
            
        density = clear_inner_disk(density)
        xs_grid, ys_grid, density_cart = polar_to_cartesian(density, rad, theta)

        # Axis
        #if axis == "zoom":
        #    x = (rad - 1) / scale_height
        #    prefix = "zoom_"
        #    plot.xlim(0, 40) # to match the ApJL paper
        #    xlabel = r"($r - r_p$) $/$ $h$"
        #else:
        #    x = rad
        #    prefix = ""
        #    plot.xlim(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]))
        #    xlabel = "Radius"
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
        result = ax.pcolormesh(xs_grid, ys_grid, np.transpose(density_cart), cmap = cmap)
        cbar = fig.colorbar(result)
        result.set_clim(clim[0], clim[1])

        cbar.set_label(r"Normalized Surface Density", fontsize = fontsize + 2, rotation = 270, labelpad = 20)

        # Get rid of interior
        circle = plot.Circle((0, 0), min(rad), color = "black")
        fig.gca().add_artist(circle)

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
        title2 = r"$t = %d$ $\rm{orbits}}$, $m_p(t) = %.2f$ $M_{Jup}$" % (orbit, current_mass)
        #plot.xlabel("x", fontsize = fontsize)
        #plot.ylabel("y", fontsize = fontsize)
        plot.title("%s" % (title2), y = 1.01, fontsize = fontsize + 1)
        plot.text(0.0, 3.14, title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

        plot.text(0.9 * plot.xlim()[0], 0.8 * plot.ylim()[1], size_label[size] + r"$\rm{-size}$", fontsize = fontsize, color = "white", bbox = dict(facecolor = 'none', edgecolor = 'white', linewidth = 1.5, pad = 7.0))

        # Save and Close
        plot.savefig("%s/%ssquareDustDensityMap_%04d_%s.png" % (save_directory, prefix, i, size), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    #choose_axis(i, "normal")
    choose_axis(i, "zoom")


##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    size = sys.argv[2]
    if frame_number == -1:
        # Plot Sample
        max_frame = 125 #util.find_max_frame()
        sample = np.linspace(0, max_frame, 126) # 100 evenly spaced frames
        #for i in sample:
        #    make_plot(i)

        p = Pool(10)
        p.map(make_plot, size, sample)
        p.terminate()

    else:
        # Plot Single
        make_plot(frame_number, size, show = True)
else:
    # Search for maximum frame
    density_files = glob.glob("gasdens*.dat")
    max_frame = find_max_frame()
    num_frames = max_frame + 1

    #for i in range(num_frames):
    #    make_plot(i)

    #### ADD TRY + CATCH BLOCK HERE!!!!! ####

    #p = Pool() # default number of processes is multiprocessing.cpu_count()
    #p.map(make_plot, range(num_frames))
    #p.terminate()

    #### Make Movies ####
    #make_movies()



