"""
plots "consecutive" (w/ rate) azithumal profiles around the peak radial density

Usage:
python plotAzimuthalDensity.py
"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool

import math
import numpy as np

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
    outer_disk_end = np.searchsorted(rad, 2.2) # look for max radial density before r = 2.6
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start:])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_density

def find_min(averagedDensity, peak_rad):
    try:
        outer_disk_start = np.searchsorted(rad, 1.0) # look for max radial density beyond r = 1.1
        outer_disk_end = np.searchsorted(rad, peak_rad)
        min_rad_outer_index = np.argmin(averagedDensity[outer_disk_start : outer_disk_end])

        min_index = outer_disk_start + min_rad_outer_index
        min_rad = rad[min_index]
        min_density = averagedDensity[min_index]

        #print "Min", min_rad, min_density
        return min_rad, min_density
    except:
        # No Gap Yet
        return peak_rad, 0

#### Data ####

vtheta_keplerian = np.array(fromfile("gasvtheta0.dat").reshape(num_rad, num_theta))

num_modes = 6
default_modes = range(1, num_modes + 1)

def get_data(frame_i, frame, modes = default_modes):
    """ frame_i is ith frame, frame is frame number """

    # Load Data Files
    normalized_density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(normalized_density, axis = 1)

    vrad = np.array(fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
    vtheta = np.array(fromfile("gasvtheta%d.dat" % frame).reshape(num_rad, num_theta)) - vtheta_keplerian

    vorticity = util.velocity_curl(vrad, vtheta, rad, theta, frame = ref_frame)
    combo = normalized_density[:-1, :-1] * vorticity

    # Find Peak in Radial Profile (in Outer Disk)
    peak_rad, peak_density = find_peak(averagedDensity)
    min_rad, min_density = find_min(averagedDensity, peak_rad)

    # Gather Azimuthal Profiles
    num_profiles = 5
    spread = 1.0 * scale_height # half-width

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = np.array([np.fft.fft(combo[azimuthal_index, :]) for azimuthal_index in azimuthal_indices])

    # Normalize by m = 0 mode (integral of density), Take Absolute Value
    azimuthal_profiles = np.array([np.abs(azimuthal_profile / azimuthal_profile[0]) for azimuthal_profile in azimuthal_profiles])

    for m, mode in enumerate(modes):
        modes_over_time[m, frame_i] = np.max(azimuthal_profiles[:, mode])

    single_mode_strength[frame_i] = modes_over_time[0, frame_i] / np.max(modes_over_time[1:, frame_i]) # m = 1 / Max of Higher Number Modes

    print "%d: %.4f, %.4f, %.4f, %.4f, %.4f - [%.4f]" % (frame, np.max(azimuthal_profiles[:, 1]), np.max(azimuthal_profiles[:, 2]), np.max(azimuthal_profiles[:, 3]), np.max(azimuthal_profiles[:, 4]), np.max(azimuthal_profiles[:, 5]), single_mode_strength[frame_i])

#### Gather Data Over Time ####

## Use These Frames ##
rate = 5 # 5 works better, but is very slow
start = 10
max_frame = util.find_max_frame()
frame_range = range(start, max_frame, rate)

## Track Modes ##
modes_over_time = np.zeros((num_modes, len(frame_range)))
single_mode_strength = np.zeros(len(frame_range))

for i, frame in enumerate(frame_range):
    get_data(i, frame)

##### PLOTTING #####

# Make Directory
#directory = "fftAzimuthalDensity"
#try:
#    os.mkdir(directory)
#except:
#    print "Directory Already Exists"

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 3

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    ### Plot ###

    for i, mode in enumerate(default_modes):
        alpha = 0.4
        if mode == 1:
            alpha = 1.0
        if mode == 3 or mode == 5:
            alpha = 0.7
        plot.plot(frame_range, modes_over_time[i, :], linewidth = linewidth, alpha = alpha, label = "%d" % mode)

    plot.plot(frame_range, single_mode_strength, color = "black", linewidth = linewidth)
    plot.plot(frame_range, np.ones(len(frame_range)), color = "black", linewidth = 1) # Reference Line at 1.0

    # Axis
    plot.xlim(0, frame_range[-1])
    plot.ylim(10**(-3.5), 10**(1.0))
    plot.yscale("log")

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel(r"$\Sigma (v - v_K)$ Mode Amplitudes", fontsize = fontsize)
    plot.title("%s" % (this_title), fontsize = fontsize + 1)

    plot.legend(loc = "upper right", bbox_to_anchor = (1.2, 1.0)) # outside of plot

    # Save and Close
    plot.savefig("fft_combination_modes.png", bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

make_plot()