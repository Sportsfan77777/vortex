"""
convolves intensity map to make it look like an alma image

python convolveIntensityMap.py intensitymap.dat
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np

from scipy import signal
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

save_directory = "intensityMaps"

# System Parameters
radius = 20.0 # radius of planet (in AU)
radius_unit = radius * (1.496 * 10**13) # (AU / cm)

# Image Parameters
beam_diameter = 10.0
beam_size = 0.5 * (beam_diameter / radius) # the sigma of the gaussian, not the beam diameter

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

rad = np.loadtxt("radial.dat") / radius_unit
theta = np.loadtxt("azimuthal.dat")
theta = np.linspace(0, 2.0 * np.pi, len(theta)) # This one is better for plotting

num_rad = len(rad)
num_theta = len(theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])
viscosity = float(fargo_par["Viscosity"])

planet_mass = float(fargo_par["PlanetMass"])
taper_time = int(float(fargo_par["MassTaper"]))

## Get RT Parameters ##
rt_par = np.loadtxt("parameters.dat")
wavelength = int(float(rt_par[2]))

### Helper Functions ###

# Converter #

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

    return xs, ys, xs_grid, ys_grid, cart_data

# Convolver #

def convolve_intensity(intensity):
    # Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.fftconvolve.html#scipy.signal.fftconvolve

    # Determine Gaussian Parameters
    dr = rad[1] - rad[0]
    sigma = int(beam_size / dr)
    window = signal.gaussian(5 * sigma, std = sigma)

    # Construct 2-D Gaussian Kernel
    normed_window = window / np.sum(window)
    kernel = np.outer(normed_window, normed_window)

    convolved_intensity = signal.fftconvolve(intensity, kernel, mode = 'same')
    return convolved_intensity

# Contraster #

def record_contrast(intensity, xs, ys):
    # Get intensity along x = 0
    zero_i = np.searchsorted(xs, 0)
    two_sliver = intensity[zero_i - 1 : zero_i + 1, :] # two columns around x = 0
    sliver = np.average(two_sliver, axis = 0)

    # Mask inner disk
    lower_i = np.searchsorted(ys, -1)
    upper_i = np.searchsorted(ys, 1)
    sliver[lower_i : upper_i] = 0.0

    # Find argmax (y-coor, and opposite y-coor)
    max_yi = np.argmax(sliver)
    max_y = ys[max_yi]

    opposite_i = np.searchsorted(ys, -max_y)

    # Mark contrast, max, opposite
    maximum = np.max(sliver)
    opposite = sliver[opposite_i]
    contrast = maximum / opposite

    return contrast, maximum, opposite

# Saver #

def save_in_polar(intensity_cart, xs, ys, order = 3):
    ### Convert to Polar ###

    # Set up rt-grid
    old_rs = rad
    old_thetas = theta

    rs_grid, thetas_grid = np.meshgrid(old_rs, old_thetas)

    # Interpolate xy-grid
    interpolated_xs = interpolate(xs, np.arange(len(xs)), bounds_error = False)
    interpolated_ys = interpolate(ys, np.arange(len(ys)), bounds_error = False)

    # Match up rt-grid with xy-grid
    new_xs = rs_grid * np.cos(thetas_grid)
    new_ys = rs_grid * np.sin(thetas_grid)

    new_interpolated_xs = interpolated_xs(new_xs.ravel())
    new_interpolated_ys = interpolated_ys(new_ys.ravel())

    # Convert (Interpolate)
    polar_data = map_coordinates(intensity_cart, np.array([new_interpolated_xs, new_interpolated_ys]), order = order).reshape(new_xs.shape)

    # Save
    save_fn = "convolved_polar_intensity.npy"
    np.save(save_fn, intensity_polar)

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "inferno"
#clim = [0, 2]
clim_in = np.searchsorted(rad, 1.1)
clim_out = np.searchsorted(rad, 2.3)

clim_low = 10
clim_high = 96

fontsize = 14
my_dpi = 100

def make_plot(frame, show = False):
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
        data = np.loadtxt("intensitymap.out")
        intensity = data[:, -1].reshape(num_rad, num_theta)
        xs, ys, xs_grid, ys_grid, intensity_cart = polar_to_cartesian(intensity, rad, theta)

        # Convolve Data
        convolved_intensity = convolve_intensity(intensity_cart)
        contrast, maximum, opposite = record_contrast(convolved_intensity, xs, ys)

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
        result = ax.pcolormesh(xs_grid, ys_grid, np.transpose(convolved_intensity), cmap = cmap)
        cbar = fig.colorbar(result)

        clim = [np.percentile(intensity[clim_in : clim_out, :], clim_low), np.percentile(intensity[clim_in : clim_out, :], clim_high)]
        result.set_clim(clim[0], clim[1])

        cbar.set_label(r"Intensity", fontsize = fontsize + 2, rotation = 270, labelpad = 20)

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

    i = frame
    #choose_axis(i, "normal")
    choose_axis(i, "zoom")


##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    if frame_number == -1:
        # Plot Sample
        max_frame = 125 #util.find_max_frame()
        sample = np.linspace(0, max_frame, 126) # 100 evenly spaced frames
        #for i in sample:
        #    make_plot(i)

        p = Pool(10)
        p.map(make_plot, sample)
        p.terminate()

    else:
        # Plot Single
        make_plot(frame_number, show = True)
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



