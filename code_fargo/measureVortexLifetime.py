"""
Measure vortex lifetime

Criteria from Fu+14:
A vortex is deemed “dead” after either the averaged azimuthal density
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

num_frames = fargo_par["Ntot"] * fargo_par["Ninterm"]

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
    return array[np.searchsorted(used_rad, start) : np.searchsorted(used_rad, stop)]

def find_density_max(frame):
    """ returns radius with maximum azimuthally averaged density """
    i = frame
    # Data
    truncated_rad = truncate(used_rad)

    density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta)) / surface_density_zero
    avg_density = truncate(np.average(density, axis = 1))

    kernel_size = len(avg_density) / 5
    smoothed_avg_density = smooth(avg_density, kernel_size)

    arg_max = np.argmax(smoothed_avg_density)
    radius_max = truncated_rad[arg_max]

    return radius_max

def find_vortensity_min(frame):
    """ returns radius with maximum azimuthally averaged density """
    pass


##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
my_dpi = 100

fontsize = 14
linewidth = 4

def plot_vortex_location(min_frame = 100, max_frame = num_frames):
    frame_range = range(min_frame, max_frame)

    vortex_locations = []
    for frame in frame_range:
        vortex_location = find_density_max(frame)
        vortex_locations.append(vortex_location)


    # Set up figure
    fig = plot.figure()
    plot.plot(frame_range, vortex_locations, linewidth = linewidth)

    # Annotate
    plot.xlabel("Orbit", fontsize = fontsize)
    plot.ylabel("Azimuthally Averaged Vortensity", fontsize = fontsize)
    plot.title("Vortex Location", fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("%s/vortexLocation_byDensity.png" % (save_directory), bbox_inches = 'tight', dpi = my_dpi)
    #plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


def make_plot(interval):
    # Plot "density variation" and "vortensity variation" around the center of the vortex


#### PLOTTING ####

# Search for maximum frame
density_files = glob.glob("gasdens*.dat")
max_frame = 0
for d_f in density_files:
    name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
    frame_number = int(name[7:]) # just 999
    if frame_number > max_frame:
        max_frame = frame_number

plot_vortex_location(max_frame = max_frame)

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


