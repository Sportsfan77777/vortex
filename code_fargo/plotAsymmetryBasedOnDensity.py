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

def gaussian(num_points, scale = 1):
    """ calculates Gaussian weights """
    def helper(x, scale = 1):
        exponent = -0.5 * x**2 / scale**2
        return np.exp(exponent)

    inputs = 1.0 * np.array(range(num_points))
    inputs -= np.average(inputs) # center on zero

    array = [helper(value) for value in inputs]
    array /= np.sum(array) # Normalize

    return array


def find_radial_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.2)
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start : outer_disk_end])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_index, peak_density

def find_azimuthal_peak(azimuthalDensity):
    twenty_degrees = (np.pi / 180.0) * 20.0
    kernel_size = np.searchsorted(rad, twenty_degrees) # kernel corresponds to 8 degrees

    smoothed_density = smooth(azimuthalDensity, kernel_size)
    peak_theta_index = np.argmax(smoothed_density)
    peak_theta = theta[peak_theta_index]

    return peak_theta, peak_theta_index

#### Data ####

threshold = 0.77 # Vortex is Above Threshold ### 1.0 / 1.3 is used because that is the initial density where the vortex forms

def measure_asymmetry(frame):
    # Find Peak in Radial Profile (in Outer Disk)
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    peak_rad, peak_rad_index, peak_density = find_radial_peak(averagedDensity)

    # Take Weighted Average of Azimuthal Profiles
    num_profiles = 7
    weights = gaussian(num_profiles, scale = num_profiles / 2.0)
    spread = 1.5 * scale_height # half-width

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = np.array([density[azimuthal_index, :] for azimuthal_index in azimuthal_indices])

    weighted_azimuthal_profile = np.average(azimuthal_profiles, weights = weights, axis = 0)

    # Find Peak in Azimuthal Profile (and center around that peak)
    peak_theta, peak_theta_index = find_azimuthal_peak(weighted_azimuthal_profile)

    initial_center = np.searchsorted(theta, np.pi)

    centered_profiles = np.roll(azimuthal_profiles, initial_center - peak_theta_index, axis = 1)
    centered_profile = np.roll(weighted_azimuthal_profile, initial_center - peak_theta_index)

    # Choose Threshold
    #pass

    # Find Bounds Marked by Threshold
    high_bound_index = np.argmax(centered_profile[initial_center : ] < threshold)
    low_bound_index = np.argmax((centered_profile[::-1])[len(centered_profile) - initial_center : ] < threshold)

    # Convert Back to Thetas and Indices
    high_theta_index = (initial_center) + high_bound_index
    low_theta_index = (len(centered_profile) - initial_center) - low_bound_index

    high_theta = theta[high_theta_index]
    low_theta = theta[low_theta_index]

    azimuthal_extent = (180.0 / np.pi) * (high_theta - low_theta)

    print "%d: %.1f, %d, %d" % (frame, azimuthal_extent, low_theta_index, high_theta_index)

    if azimuthal_extent > 5:
        # Find Quartiles (25%, 50%, 75%)
        lower_quartile = np.percentile(centered_profiles[:, low_theta_index : high_theta_index], 25)
        avg_density = np.percentile(centered_profiles[:, low_theta_index : high_theta_index], 50)
        upper_quartile = np.percentile(centered_profiles[:, low_theta_index : high_theta_index], 75)

        return azimuthal_extent, avg_density, lower_quartile, upper_quartile
    else:
        # No Vortex Yet (or any features really)
        return azimuthal_extent, threshold, threshold, threshold

## Use These Frames ##
rate = 25 # 5 works better, but is very slow
start_of_vortex = 700
max_frame = 1500 #util.find_max_frame()
frame_range = range(start_of_vortex, max_frame, rate)

vortex_azimuthal_widths = []

vortex_avg_densities = []
upper_quartiles = []
lower_quartiles = []

for frame in frame_range:
    asymmetry, avg_density, lower_quartile, upper_quartile = measure_asymmetry(frame)

    vortex_azimuthal_widths.append(asymmetry)

    vortex_avg_densities.append(avg_density)
    lower_quartiles.append(lower_quartile)
    upper_quartiles.append(upper_quartile)


## Smooth Each Array
kernel_size = 5

smoothed_vortex_azimuthal_widths = smooth(vortex_azimuthal_widths, kernel_size)
smoothed_vortex_avg_densities = smooth(vortex_avg_densities, kernel_size)

smoothed_upper_quartiles = smooth(upper_quartiles, kernel_size)
smoothed_lower_quartiles = smooth(lower_quartiles, kernel_size)

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

color = ["blue", "red", "green"]

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ### Hatched Region for Quartiles ###
    ax2.fill_between(frame_range, smoothed_lower_quartiles, smoothed_upper_quartiles, color = color[1], facecolor = "None", alpha = normal_alpha * offset / 1.3, hatch = "\\")

    ### Plot ###
    ax2.plot(frame_range, vortex_avg_densities, color = color[1], linewidth = linewidth - 2, alpha = normal_alpha * offset)
    ax2.plot(frame_range, smoothed_vortex_avg_densities, color = color[1], linewidth = linewidth - 1, alpha = smooth_alpha * offset)

    ax1.plot(frame_range, vortex_azimuthal_widths, color = color[0], linewidth = linewidth - 2, alpha = normal_alpha)
    ax1.plot(frame_range, smoothed_vortex_azimuthal_widths, color = color[0], linewidth = linewidth, alpha = smooth_alpha) # Dominant Line (that is why it is last)

    # Limits
    plot.xlim(frame_range[0], frame_range[-1])

    #if np.max(smoothed_vortex_azimuthal_widths) > 180:
    #    angles = np.linspace(0, 360, 7)
    #    ax1.set_ylim(0, 360)
    #    ax1.set_yticks(angles)
    #else:
    #    angles = np.linspace(0, 180, 7)
    #    ax1.set_ylim(0, 180)
    #    ax1.set_yticks(angles)
    angles = np.linspace(0, 360, 7)
    ax1.set_ylim(0, 360)
    ax1.set_yticks(angles)

    max_density = np.max(smoothed_upper_quartiles)
    max_y = np.ceil(2.0 * max_density) / 2.0 # round up to the nearest 0.5
    ax2.set_ylim(threshold, max_y)

    # Annotate
    this_title = readTitle()
    ax1.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
    ax1.set_ylabel("Azimuthal Extent", fontsize = fontsize, color = color[0])
    ax2.set_ylabel("Mean Density", fontsize = fontsize, color = color[1])
    plot.title("Vortex Properties: %s" % (this_title), fontsize = fontsize + 1)

    for tick1 in ax1.get_yticklabels():
        tick1.set_color(color[0])

    for tick2 in ax2.get_yticklabels():
        tick2.set_color(color[1])

    # Save and Close
    plot.savefig("vortex_asymmetry77_v2.png", bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


# Plot!
make_plot()

