"""
plots asymmtry and mean density of vortex around the peak radial density over time

Usage:
python plotAsymmetry.py
"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool

import math
import numpy as np
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import rcParams # replace with rc ???
from pylab import fromfile

import util
from readTitle import readTitle

## Check frame ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    ref_frame = 0
else:
    # fargo
    ref_frame = 1

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

surface_density = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

### Helper Methods ###
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

def find_radial_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.2)
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start : outer_disk_end])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_index, peak_density

def find_azimuthal_extrema(azimuthalProfile, maximum = True):
    ### Use Smoothed Profile ###
    if maximum:
        extrema = np.argmax
    else:
        extrema = np.argmin

    twenty_degrees = (np.pi / 180.0) * 20.0
    kernel_size = np.searchsorted(rad, twenty_degrees) # kernel corresponds to 20 degrees

    smoothed_profile = smooth(azimuthalProfile, kernel_size)
    peak_theta_index = extrema(smoothed_profile)

    peak_theta = theta[peak_theta_index]
    peak_value = azimuthalProfile[peak_theta_index]

    return peak_theta, peak_theta_index, peak_value

#### Data ####

def find_extrema(i, frame):
    # Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    vrad = (fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % frame).reshape(num_rad, num_theta))

    vorticity = util.velocity_curl(vrad, vtheta, rad, theta, frame = ref_frame)
    vortensity = vorticity / density[1:, 1:]

    # Find Peak in Radial Profile (in Outer Disk)
    peak_rad, peak_rad_index, peak_density = find_radial_peak(averagedDensity)

    # Find Peaks in Azimuthal Profiles (and center around that peak)
    density_peak_theta, density_theta_index, max_density = find_azimuthal_extrema(density[peak_rad_index], maximum = True) # Max
    vortensity_peak_theta, vortensity_peak_theta_index, min_vortensity = find_azimuthal_extrema(vortensity[peak_rad_index], maximum = False) # Min

    print "%d, %.2f, %.3f" % (frame, max_density, min_vortensity)

    #return max_density, min_vortensity
    maximum_densities[i] = max_density
    minimum_vortensities[i] = min_vortensity

## Use These Frames ##
rate = 5 # 5 works better, but is very slow
start_of_vortex = 0
max_frame = util.find_max_frame()
frame_range = range(start_of_vortex, max_frame, rate)

maximum_densities = np.array(len(frame_range))
minimum_vortensities = np.array(len(frame_range))

for i, frame in frame_range:
    find_extrema(i, frame)

## Smooth Each Array ##
kernel_size = 5

smoothed_maximum_densities = smooth(maximum_densities, kernel_size)
smoothed_minimum_vortensities = smooth(minimum_vortensities, kernel_size)

## Store Max Density ##
peak_density = np.max(smoothed_maximum_densities)
pickle.dump(peak_density, open("peak_density.p"))

##### PLOTTING #####

# Make Directory
#directory = "asymmetry"
#try:
#    os.mkdir(directory)
#except:
#    print "Directory Already Exists"

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

smooth_alpha = 1.0
normal_alpha = 0.35
offset = 0.7 # 70%

fontsize = 14
linewidth = 4

color = ["blue", "red"]

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ### Plot ###
    ax2.plot(frame_range, minimum_vortensities, color = color[1], linewidth = linewidth - 2, alpha = normal_alpha * offset)
    ax2.plot(frame_range, smoothed_minimum_vortensities, color = color[1], linewidth = linewidth - 1, alpha = smooth_alpha * offset)

    ax1.plot(frame_range, maximum_densities, color = color[0], linewidth = linewidth - 2, alpha = normal_alpha)
    ax1.plot(frame_range, smoothed_maximum_densities, color = color[0], linewidth = linewidth, alpha = smooth_alpha) # Dominant Line (that is why it is last)

    # Limits
    plot.xlim(frame_range[0], frame_range[-1])

    ax1.set_ylim(0, 3.5) # Density Range
    ax2.set_ylim(0, 0.3) # Vortensity Range

    # Annotate
    this_title = readTitle()
    ax1.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
    ax1.set_ylabel("Maximum Density", fontsize = fontsize, color = color[0])
    ax2.set_ylabel("Minimum Vortensity", fontsize = fontsize, color = color[1])
    plot.title("Vortex Properties: %s" % (this_title), fontsize = fontsize + 1)

    for tick1 in ax1.get_yticklabels():
        tick1.set_color(color[0])

    for tick2 in ax2.get_yticklabels():
        tick2.set_color(color[1])

    # Save and Close
    plot.savefig("vortex_extrema.png", bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

# Plot!
make_plot()

