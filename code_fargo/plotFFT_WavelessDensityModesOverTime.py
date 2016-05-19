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
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot
from matplotlib import gridspec

from itertools import groupby

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

def find_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.2) # look for max radial density before r = 2.6
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start : outer_disk_end])

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

def mark_vortex_start(frame_range, single_mode, single_mode_strength):
    # Find first frame beyond 50 (not too early) where single_mode > 0.1 (vortex is strong) and single_mode_strength > 1.0 (vortex is dominant)
    start = -1
    for (frame, single_mode_i, single_mode_strength_i) in zip(frame_range, single_mode, single_mode_strength):
        if frame > 50 and single_mode_i > 0.1 and single_mode_strength_i > 1.0:
            start = frame
            break

    print "Vortex Start: %d" % start

    # Store in a pickle
    start_fn = "start.p"
    pickle.dump(start, open(start_fn, "wb"))

def mark_vortex_end(frame_range, single_mode, single_mode_strength, single_mode_concentration):
    # Helper
    def find_consecutive_ranges(array, values, cutoff, ranges = [], greater = True):
        if greater:
            test = lambda x : values[np.searchsorted(array, x)] > cutoff
        else:
            test = lambda x : values[np.searchsorted(array, x)] < cutoff

        for (match, group) in groupby(array, key = test):
            # Identify start and end of each range
            if match:
                start = next(group)
                end = start
                for end in group:
                    pass
                this_range = [start, end]
                ranges.append(this_range)

        return ranges

    cutoff = 0.1
    ranges = find_consecutive_ranges(frame_range, single_mode, cutoff)

    cutoff = 0.1
    ranges = find_consecutive_ranges(frame_range, single_mode_concentration, cutoff, ranges = ranges, greater = False)

    cutoff = 1.0
    ranges = find_consecutive_ranges(frame_range, single_mode_strength, cutoff, ranges = ranges)

    print "Vortex End Candidates: ", [r[-1] for r in ranges]

    # Store in a pickle
    end_candidates_fn = "end_candidates.p"
    pickle.dump(ranges, open(end_candidates_fn, "wb"))

def subtract_wave(density):
    averagedDensity = np.average(density, axis = 1)
    peak_rad, peak_density = find_peak(averagedDensity)

    end_radius = min(peak_rad + 7 * scale_height, 2.5)

    used_rad, wave_locations = util.getWaveLocation(density, rad, theta, end_radius = end_radius)

    # Interpolate in Wave Region
    delta_theta = 17.5 * (np.pi / 180.0)

    for i, (r_i, wave_location) in enumerate(zip(used_rad, wave_locations)):
        rad_index = np.searchsorted(rad, r_i)

        lower_theta = wave_location - delta_theta
        upper_theta = wave_location + delta_theta

        if lower_theta < 0.0:
            lower_theta_i = np.searchsorted(theta, lower_theta + 2 * np.pi)
            upper_theta_i = np.searchsorted(theta, upper_theta)

            interpolated_thetas = np.concatenate((theta[lower_theta_i :] - 2 * np.pi, theta[ : upper_theta_i]))
            thetas = np.array([lower_theta, upper_theta])
            values = np.array([density[rad_index, lower_theta_i], density[rad_index, upper_theta_i]])

            interpolated_values = np.interp(interpolated_thetas, thetas, values)
            density[rad_index, lower_theta_i :] = interpolated_values[: len(theta[lower_theta_i :])]
            density[rad_index, : upper_theta_i] = interpolated_values[len(theta[len_theta_i :]) :]

        elif upper_theta > 2 * np.pi:
            lower_theta_i = np.searchsorted(theta, lower_theta)
            upper_theta_i = np.searchsorted(theta, upper_theta - 2 * np.pi)

            interpolated_thetas = np.concatenate((theta[lower_theta_i :], theta[ : upper_theta_i] + 2 * np.pi))
            thetas = np.array([lower_theta, upper_theta])
            values = np.array([density[rad_index, lower_theta_i], density[rad_index, upper_theta_i]])

            interpolated_values = np.interp(interpolated_thetas, thetas, values)
            density[rad_index, lower_theta_i :] = interpolated_values[: len(theta[lower_theta_i :])]
            density[rad_index, : upper_theta_i] = interpolated_values[len(theta[lower_theta_i :]) :]

        else:
            lower_theta_i = np.searchsorted(theta, lower_theta)
            upper_theta_i = np.searchsorted(theta, upper_theta)

            interpolated_thetas = theta[lower_theta_i : upper_theta_i]
            thetas = np.array([lower_theta_i, upper_theta_i])
            values = np.array([density[rad_index, lower_theta_i], density[rad_index, upper_theta_i]])

            density[rad_index, lower_theta_i : upper_theta_i] = np.interp(interpolated_thetas, thetas, values)

    return density

#### Data ####

num_modes = 6
default_modes = range(1, num_modes + 1)

def get_data(frame_i, frame, modes = default_modes):
    """ frame_i is ith frame, frame is frame number """

    # Find Peak in Radial Profile (in Outer Disk)
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    density = subtract_wave(density) # First, subtract wave!

    peak_rad, peak_density = find_peak(averagedDensity)
    min_rad, min_density = find_min(averagedDensity, peak_rad)

    # Gather Azimuthal Profiles
    num_profiles = 7
    spread = 1.0 * scale_height # half-width

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = np.array([np.fft.fft(density[azimuthal_index, :]) for azimuthal_index in azimuthal_indices])

    # Normalize by m = 0 mode (integral of density), Take Absolute Value
    azimuthal_profiles = np.array([np.abs(azimuthal_profile / azimuthal_profile[0]) for azimuthal_profile in azimuthal_profiles])

    for m, mode in enumerate(modes):
        modes_over_time[m, frame_i] = np.max(azimuthal_profiles[:, mode]) #np.sqrt(np.sum(np.power(azimuthal_profiles[:, mode], 2))) 

    single_mode_strength[frame_i] = modes_over_time[0, frame_i] / np.max(modes_over_time[1:, frame_i]) # m = 1 / Max of Higher Number Modes
    single_mode_concentration[frame_i] = np.std(azimuthal_profiles[:, 1]) / modes_over_time[0, frame_i]

    print "%d: %.4f, %.4f, %.4f, %.4f, %.4f - [%.4f, %.4f]" % (frame, np.max(azimuthal_profiles[:, 1]), np.max(azimuthal_profiles[:, 2]), np.max(azimuthal_profiles[:, 3]), np.max(azimuthal_profiles[:, 4]), np.max(azimuthal_profiles[:, 5]), single_mode_strength[frame_i], single_mode_concentration[frame_i])

#### Gather Data Over Time ####

## Use These Frames ##
rate = 5 # 5 works better, but is very slow
start = 10
max_frame = util.find_max_frame()
frame_range = np.array(range(start, max_frame + 1, rate))

## Track Modes ##
modes_over_time = np.zeros((num_modes, len(frame_range)))
single_mode_strength = np.zeros(len(frame_range))
single_mode_concentration = np.zeros(len(frame_range))

for i, frame in enumerate(frame_range):
    get_data(i, frame)

# Smooth?
kernel_size = 5
#single_mode_strength = smooth(single_mode_strength, kernel_size)
single_mode_concentration = smooth(single_mode_concentration, kernel_size)

# Highlight Vortex

vortex_highlighter = np.copy(single_mode_strength)
vortex_highlighter[vortex_highlighter < 1] -= 10**5 # m = 1 subdominant, make negative to remove from log plot
#vortex_highlighter[modes_over_time[0, :] < 0.1] -= 10**5 # m = 1 < 0.1

# Mark Start and End
single_mode = modes_over_time[0]

mark_vortex_start(frame_range, single_mode, single_mode_strength)
mark_vortex_end(frame_range, single_mode, single_mode_strength, single_mode_concentration)

##### PLOTTING #####

# Make Directory
directory = "fftAzimuthalDensity"
try:
    os.mkdir(directory)
except:
    print "Directory Already Exists"

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 3

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    gs = gridspec.GridSpec(2, 1, height_ratios = [5, 11])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex = ax1)

    #fig.subplots_adjust(hspace = 1)

    ### Plot ###

    for i, mode in enumerate(default_modes):
        alpha = 0.4
        if mode == 1:
            alpha = 1.0
        if mode == 3 or mode == 5:
            alpha = 0.7
        ax2.plot(frame_range, modes_over_time[i, :], linewidth = linewidth, alpha = alpha, label = "%d" % mode)

    ax2.plot([0, frame_range[-1]], 0.1 * np.ones(2), color = "black", linewidth = 1) # Reference Line at 0.1

    ax1.plot(frame_range, single_mode_strength, color = "black", linewidth = linewidth)
    #ax1.plot(frame_range, vortex_highlighter, color = "orange", linewidth = linewidth)
    ax1.plot([0, frame_range[-1]], np.ones(2), color = "black", linewidth = 1) # Reference Line at 1.0

    ax1.plot(frame_range, single_mode_concentration, color = "red", linewidth = linewidth - 1)
    ax1.plot([0, frame_range[-1]], 0.1 * np.ones(2), color = "black", linewidth = 1) # Reference Line at 0.1

    # Limits
    ax2.set_xlim(0, frame_range[-1])
    ax2.set_ylim(10**(-3.5), 10**(0.0))
    ax2.set_yscale("log")

    ax1.set_ylim(10**(-1.5), 10**(1.0))
    ax1.set_yscale("log")

    # Annotate
    this_title = readTitle()
    ax2.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
    ax2.set_ylabel("Density Mode Amplitudes", fontsize = fontsize)
    ax1.set_ylabel("m = 1 Strength", fontsize = fontsize)
    ax1.set_title("%s" % (this_title), fontsize = fontsize + 1)

    ax2.legend(loc = "upper right", bbox_to_anchor = (1.2, 1.0)) # outside of plot

    # Save and Close
    plot.savefig("fft_waveless_density_modes.png", bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

make_plot()