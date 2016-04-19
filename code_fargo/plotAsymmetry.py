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
def find_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start:])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_density

#### Data ####

threshold = 1.0 # Vortex is Above Threshold

def measure_asymmetry(frame):
    # Find Peak in Radial Profile (in Outer Disk)
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    peak_rad, peak_density = find_peak(averagedDensity)

    # Gather Azimuthal Profiles
    num_profiles = 5
    spread = 1.0 * scale_height # half-width

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = [density[azimuthal_index, :] for azimuthal_index in azimuthal_indices]

    # For each profile, measure the azimuthal extent of the vortex (density > threshold = 1.0)
    vortex_value = 360.0 / float(fargo_par["Nsec"])
    background_value = 0.0

    asymmetry_values = []
    avg_densities = []
    for azimuthal_profile in azimuthal_profiles:
        # Get Asymmetry
        copy = np.zeros(np.shape(azimuthal_profile))

        copy[azimuthal_profile < threshold] = background_value
        copy[azimuthal_profile >= threshold] = vortex_value

        azithumal_extent = np.sum(copy) # should be the angle spanned by the vortex (in degrees)
        asymmetry_values.append(azithumal_extent)

        # Get Density
        copy = azimuthal_profile[azimuthal_profile > threshold]

        avg_density = np.mean(copy)
        avg_densities.append(avg_density)

    # Take Median of Profiles between -H and +H about peak in radial density profile
    asymmetry = np.median(asymmetry_values)

    # Get Mean Density of Vortex within One Scale Height
    start = azimuthal_indices[0]
    end = azimuthal_indices[-1]
    vortex_zone = density[start : end, :]

    vortex_densities = vortex_zone[vortex_zone > threshold]
    avg_density = np.mean(vortex_densities)

    return asymmetry, avg_density

## Use These Frames ##
rate = 10
start_of_vortex = 10
max_frame = util.find_max_frame()
frame_range = range(start_of_vortex, max_frame, rate)

vortex_azimuthal_widths = []
vortex_avg_densities = []

for frame in frame_range:
    asymmetry, avg_density = measure_asymmetry(frame)

    vortex_azimuthal_widths.append(asymmetry)
    vortex_avg_densities.append(avg_density)

## Smooth Both Arrays
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

kernel_size = 10
vortex_azimuthal_widths = smooth(vortex_azimuthal_widths, kernel_size)
vortex_avg_densities = smooth(vortex_avg_densities, kernel_size)

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

alpha = 0.65
fontsize = 14
linewidth = 4

color = ["blue", "red"]

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ### Plot ###
    ax2.plot(frame_range, vortex_avg_densities, color = color[1], linewidth = linewidth - 1, alpha = alpha)
    ax1.plot(frame_range, vortex_azimuthal_widths, color = color[0], linewidth = linewidth) # Dominant Line (that is why it is last)

    # Limits
    if np.max(vortex_azimuthal_widths) > 180:
        ax1.set_ylim(0, 360)
    else:
        ax1.set_ylim(0, 180)

    max_density = np.max(vortex_avg_densities)
    max_y = np.ceil(2.0 * max_density) / 2.0 # round up to the nearest 0.5
    ax2.set_ylim(threshold, max_y)

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    ax1.set_ylabel("Azimuthal Extent", fontsize = fontsize, color = color[0])
    ax2.set_ylabel("Mean Density", fontsize = fontsize, color = color[1])
    plot.title("Vortex Properties: %s" % (this_title), fontsize = fontsize + 1)

    for tick1 in ax1.get_yticklabels():
        tick1.set_color(color[0])

    for tick2 in ax2.get_yticklabels():
        tick2.set_color(color[1])

    # Save and Close
    plot.savefig("vortex_asymmetry.png", bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


# Plot!
make_plot()

