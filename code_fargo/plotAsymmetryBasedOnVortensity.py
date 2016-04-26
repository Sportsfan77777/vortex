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
def find_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.2)
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start : outer_disk_end])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_density

#### Data ####

threshold = 0.25 # Vortex is Below Threshold

def measure_asymmetry(frame):
    print frame

    # Load Data Files
    normalized_density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(normalized_density, axis = 1)

    vrad = (fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
    vtheta = (fromfile("gasvtheta%d.dat" % frame).reshape(num_rad, num_theta))

    vorticity = util.velocity_curl(vrad, vtheta, rad, theta, frame = ref_frame)
    vortensity = vorticity / normalized_density[1:, 1:]

    # Find Peak in Radial Profile (in Outer Disk)

    peak_rad, peak_density = find_peak(averagedDensity)

    # Gather Azimuthal Profiles
    num_profiles = 5
    spread = 1.0 * scale_height # half-width

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = [vortensity[azimuthal_index, :] for azimuthal_index in azimuthal_indices]

    # For each profile, measure the azimuthal extent of the vortex (vortensity < threshold = 0.25)
    vortex_value = 360.0 / float(fargo_par["Nsec"])
    background_value = 0.0

    asymmetry_values = []
    avg_densities = []
    for azimuthal_profile in azimuthal_profiles:
        # Get Asymmetry
        copy = np.zeros(np.shape(azimuthal_profile))

        copy[azimuthal_profile > threshold] = background_value
        copy[azimuthal_profile <= threshold] = vortex_value

        azithumal_extent = np.sum(copy) # should be the angle spanned by the vortex (in degrees)
        asymmetry_values.append(azithumal_extent)

        # Get Vortensity
        copy = azimuthal_profile[azimuthal_profile <= threshold]

        avg_vortensity = np.mean(copy)
        avg_densities.append(avg_vortensity)

    # Take Median of Profiles between -H and +H about peak in radial density profile
    asymmetry = np.median(asymmetry_values)

    # Get Mean Density of Vortex within One Scale Height (+ 25% and 75% values)
    start = azimuthal_indices[0]
    end = azimuthal_indices[-1]
    vortex_zone = vortensity[start : end, :]

    vortex_vortensities = vortex_zone[vortex_zone < threshold]
    avg_vortensity = np.mean(vortex_vortensities)

    # Detect Nan
    if avg_vortensity != avg_vortensity:
        avg_vortensity = threshold
        lower_quartile = threshold
        upper_quartile = threshold
    else:
        # If no nan, compute quartiles
        lower_quartile = np.percentile(vortex_vortensities, 25) # 25%
        upper_quartile = np.percentile(vortex_vortensities, 75) # 75%

    return asymmetry, avg_vortensity, lower_quartile, upper_quartile

## Use These Frames ##
rate = 25 # 5 works better, but is very slow
start_of_vortex = 10
max_frame = util.find_max_frame()
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
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

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
    ax3 = ax1.twinx()

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

    ax2.set_ylim(0.0, threshold)

    # Annotate
    this_title = readTitle()
    ax1.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
    ax1.set_ylabel("Azimuthal Extent", fontsize = fontsize, color = color[0])
    ax2.set_ylabel("Mean Vortensity", fontsize = fontsize, color = color[1])
    plot.title("Vortex Properties: %s" % (this_title), fontsize = fontsize + 1)

    for tick1 in ax1.get_yticklabels():
        tick1.set_color(color[0])

    for tick2 in ax2.get_yticklabels():
        tick2.set_color(color[1])

    # Save and Close
    plot.savefig("vortex_vortensity_asymmetry.png", bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


# Plot!
make_plot()

