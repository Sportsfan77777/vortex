"""
Measure vortex lifetime

Criteria from Fu+14:
A vortex is deemed "dead" after either the averaged azimuthal density
variation or the averaged azimuthal potential vorticity variation within 10H
(scale height) wide band around the vortex drops below 10%.

Method:
(1) Find vortex center --- use density max and vortensity min as a guide
(2) 

python measureVortexLifetime.py
python measureVortexLifetime.py interval
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np

from scipy import signal as sig
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile


save_directory = "lifetime"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

num_frames = int(fargo_par["Ntot"]) * int(fargo_par["Ninterm"])

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

# Curl function
def curl(v_rad, v_theta, rad, theta):
    """ z-component of the curl (because this is a 2-D simulation)"""
    ### Start Differentials ###
    d_rad = np.diff(rad)
    d_theta = np.diff(theta)

    dv_rad = np.diff(v_rad, axis = 1)
    dv_theta = np.diff(rad[:, None] * v_theta, axis = 0)
    ### End Differentials ###

    # z-Determinant
    partial_one = dv_theta / d_rad[:, None]
    partial_two = dv_rad / d_theta

    z_curl = (partial_one[:, 1:] - partial_two[1:, :]) / rad[1:, None]

    # Shift out of rotating frame (http://arxiv.org/pdf/astro-ph/0605237v2.pdf)
    z_curl += 2

    return z_curl

### Find Vortex ###

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
# Truncate
def truncate(array, start = 1.2, stop = 3.0):
    """ truncates azimuthally averaged array between two radii """
    return array[np.searchsorted(rad, start) : np.searchsorted(rad, stop)]

def find_density_max(frame):
    """ returns radius with maximum azimuthally averaged density """
    i = frame
    # Data
    truncated_rad = truncate(rad)

    density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta)) / surface_density_zero
    avg_density = truncate(np.average(density, axis = 1))

    kernel_size = 1
    smoothed_avg_density = smooth(avg_density, kernel_size)

    # Retrieve radius with max value
    arg_max = np.argmax(smoothed_avg_density)
    radius_max = truncated_rad[arg_max]

    return radius_max

def find_vortensity_min(frame):
    """ returns radius with maximum azimuthally averaged density """
    i = frame

    # Find density max (Search for vortensity minimum only near the density max)
    density_max = find_density_max(frame)
    start = density_max - (0.2 * (scale_height / 0.06))
    stop = density_max + (0.1 * (scale_height / 0.06))

    # Data
    truncated_rad = truncate(rad, start = start, stop = stop)

    density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    vrad = (fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta))

    vorticity = curl(vrad, vtheta, rad, theta)
    vortensity = vorticity / normalized_density[1:, 1:]
    avg_vortensity = truncate(np.average(vortensity, axis = 1), start = start, stop = stop)

    kernel_size = 1
    smoothed_avg_vortensity = smooth(avg_vortensity, kernel_size)

    # Retrieve radius with min value
    arg_min = np.argmin(smoothed_avg_vortensity)
    radius_min = truncated_rad[arg_min]

    return radius_min


def get_density_variation(frame, vortex_location, figure = False):
    """
    return variation in density (defined to be max / avg at a particular radius)
    """
    i = frame
    # Data
    density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    # Azimuthal Profile
    arg_vortex = np.searchsorted(rad, vortex_location)
    kernel_size = 1 # num_theta / 200
    density_at_vortex = smooth(normalized_density[arg_vortex,:], kernel_size)

    # Plot
    if figure:
        fig = plot.figure()

        # Axis
        angles = np.linspace(0, 2 * np.pi, 7)
        degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

        plot.xlim(0, 2 * np.pi)
        plot.xticks(angles, degree_angles)

        # Plot
        plot.plot(theta, normalized_density[arg_vortex,:], linewidth = linewidth)

        # Annotate
        plot.xlabel("Theta", fontsize = fontsize)
        plot.ylabel("Density", fontsize = fontsize)

        # Save and Close
        directory = "azimuthalDensity"
        plot.savefig("%s/azimuthalDensity_%04d.png" % (directory, i), bbox_inches = 'tight', dpi = my_dpi)
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    # Variation
    max_density = np.max(density_at_vortex)
    avg_density = np.average(density_at_vortex)

    variation = max_density / avg_density
    return variation

def get_vortensity_variation(frame, vortex_location, figure = False):
    """
    return variation in vortensity (defined to be min / avg at a particular radius)
    """
    i = frame
    # Data
    density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
    normalized_density = density / surface_density_zero

    vrad = (fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta))

    vorticity = curl(vrad, vtheta, rad, theta)
    vortensity = vorticity / normalized_density[1:, 1:]

    # Azimuthal Profile
    arg_vortex = np.searchsorted(rad, vortex_location)
    kernel_size = num_theta / 200
    vortensity_at_vortex = smooth(vortensity[arg_vortex,:], kernel_size)

    # Plot
    if figure:
        fig = plot.figure()

        # Axis
        angles = np.linspace(0, 2 * np.pi, 7)
        degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

        plot.xlim(0, 2 * np.pi)
        plot.xticks(angles, degree_angles)

        # Plot
        plot.plot(theta[1:], vortensity[arg_vortex,:], linewidth = linewidth)

        # Annotate
        plot.xlabel("Theta", fontsize = fontsize)
        plot.ylabel("Vortensity", fontsize = fontsize)

        # Save and Close
        directory = "azimuthalVortensity"

        plot.savefig("%s/azimuthalVortensity_%04d.png" % (directory, i), bbox_inches = 'tight', dpi = my_dpi)
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    # Variation
    min_vortensity = np.min(vortensity_at_vortex)
    avg_vortensity = np.average(vortensity_at_vortex)

    #print min_vortensity

    variation = (1 + avg_vortensity - min_vortensity) / avg_vortensity
    return variation

##### PLOTTING #####

# Make Directories
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

try:
    os.mkdir("azimuthalDensity")
except:
    pass

try:
    os.mkdir("azimuthalVortensity")
except:
    pass


# Plot Parameters
my_dpi = 100

fontsize = 14
linewidth = 4

def plot_vortex_location(min_frame = 100, max_frame = num_frames, rate = 5, figure = True):
    """
    boolean option to plot vortex location
    regardless, return vortex location using local min of vortensity
    """

    frame_range = range(min_frame, max_frame, rate)

    vortex_locations_d = []
    vortex_locations_v = []
    for frame in frame_range:
        vortex_location_d = find_density_max(frame)
        vortex_locations_d.append(vortex_location_d)

        vortex_location_v = find_vortensity_min(frame)
        vortex_locations_v.append(vortex_location_v)

    if figure:
        # Set up figure
        fig = plot.figure()
        plot.plot(frame_range, vortex_locations_d, "b", linewidth = linewidth)
        plot.plot(frame_range, vortex_locations_v, "r", linewidth = linewidth)

        # Annotate
        plot.xlabel("Orbit", fontsize = fontsize)
        plot.ylabel("Vortex Location", fontsize = fontsize)
        #plot.title("Vortex Location", fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/vortexLocation.png" % (save_directory), bbox_inches = 'tight', dpi = my_dpi)
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    return vortex_locations_v


def plot_azimuthal_variation(vortex_locations, min_frame = 100, max_frame = num_frames, rate = 5, figure = True):
    """ 
    Plot "density variation" and "vortensity variation" around the center of the vortex
    Requires same frame range as with vortex locations
    """

    frame_range = range(min_frame, max_frame, rate)

    variations_d = []
    variations_v = []
    for i, frame in enumerate(frame_range):
        variation_d = get_density_variation(frame, vortex_locations[i])
        variations_d.append(variation_d)

        variation_v = get_vortensity_variation(frame, vortex_locations[i])
        variations_v.append(variation_v)

    # Set up figure
    fig = plot.figure()
    plot.plot(frame_range, variations_d, "b", linewidth = linewidth)
    plot.plot(frame_range, variations_v, "r", linewidth = linewidth)
    #plot.scatter(frame_range, vortex_locations)

    # Annotate
    plot.xlabel("Orbit", fontsize = fontsize)
    plot.ylabel("Variation", fontsize = fontsize)
    #plot.title("Vortex Location", fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("%s/vortexStrength.png" % (save_directory), bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

#### PLOTTING ####

# Search for maximum frame
density_files = glob.glob("gasdens*.dat")
max_frame = 0
for d_f in density_files:
    name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
    frame_number = int(name[7:]) # just 999
    if frame_number > max_frame:
        max_frame = frame_number

vortex_locations = plot_vortex_location(max_frame = max_frame)
plot_azimuthal_variation(vortex_locations, max_frame = max_frame)

# if len(sys.argv) > 1:
#     frame_number = int(sys.argv[1])
#     make_plot(frame_number)
# else:
#     # Search for maximum frame
#     density_files = glob.glob("gasdens*.dat")
#     max_frame = 0
#     for d_f in density_files:
#         name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
#         frame_number = int(name[7:]) # just 999
#         if frame_number > max_frame:
#             max_frame = frame_number
#     num_frames = max_frame + 1

#     #for i in range(num_frames):
#     #    make_plot(i)

#     p = Pool() # default number of processes is multiprocessing.cpu_count()
#     p.map(make_plot, range(num_frames))
#     p.terminate()

#     #### Make Movies ####
#     make_movies()


