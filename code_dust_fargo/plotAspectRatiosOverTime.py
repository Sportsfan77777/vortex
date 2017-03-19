"""
plots quartiles (e.g. 50%, 75%, 90%, 99%) of density in vicinity of vortex over time

Usage:
python plotQuartilesOverTime.py
"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool
from multiprocessing import Array as mp_array

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

rad = np.loadtxt("used_rad.dat")[:-1, 0]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

taper_time = int(float(fargo_par["MassTaper"]))

### Helper Functions ###

def find_peak(averagedDensity):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max radial density before r = 2.3
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start : outer_disk_end])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_density

### Data ###

def get_extents(args):
    # Unwrap Args
    i, frame = args

    # Print Update
    print "%d" % (frame)

    # Get Dust Data
    density = (fromfile("gasddens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero

    ####### Azimuthal Extents #######

    # Get Vortex Vicinity Indices
    averagedDensity = np.average(density, axis = 1)
    peak_rad, peak_density = find_peak(averagedDensity)

    vortex_start = np.max([1.05, peak_rad - 0.5 * scale_height])
    vortex_end = peak_rad + 0.5 * scale_height

    vortex_start_i = np.searchsorted(rad, vortex_start)
    vortex_end_i = np.searchsorted(rad, vortex_end)

    # Extract Near Vortex
    vortex_rad = rad[vortex_start_i : vortex_end_i]
    vortex_density = density[vortex_start_i : vortex_end_i]

    # Count number of cells above each threshold
    for t_i, threshold in enumerate(thresholds):
        # Store Copy
        tmp_vortex_density = np.copy(vortex_density)

        # Mark Above/Below Threshold
        tmp_vortex_density[tmp_vortex_density <= threshold] = 0
        tmp_vortex_density[tmp_vortex_density >  threshold] = 1

        counts = np.sum(tmp_vortex_density, axis = 1)
        median_count = np.median(counts)

        # Store Median (converted to angle)
        median_count *= (360.0 / num_theta)

        if t_i == 0:
            az_widths0[i] = median_count
        elif t_i == 1:
            az_widths1[i] = median_count
        elif t_i == 2:
            az_widths2[i] = median_count
        elif t_i == 3:
            az_widths5[i] = median_count
        else:
            az_widths10[i] = median_count

    ####### Radial Extents #######

    # Get Vortex Vicinity Indices
    averagedDensity = np.average(density, axis = 1)
    peak_rad, peak_density = find_peak(averagedDensity)

    vortex_start = np.max([1.05, peak_rad - 5.0 * scale_height])
    vortex_end = peak_rad + 5.0 * scale_height

    vortex_start_i = np.searchsorted(rad, vortex_start)
    vortex_end_i = np.searchsorted(rad, vortex_end)

    # Extract Near Vortex
    vortex_rad = rad[vortex_start_i : vortex_end_i]
    vortex_density = density[vortex_start_i : vortex_end_i]

    # Count number of cells above each threshold
    for t_i, threshold in enumerate(thresholds):
        # Store Copy
        tmp_vortex_density = np.copy(vortex_density)

        # Mark Above/Below Threshold
        tmp_vortex_density[tmp_vortex_density <= threshold] = 0
        tmp_vortex_density[tmp_vortex_density >  threshold] = 1

        counts = np.sum(tmp_vortex_density, axis = 0)
        high_count = np.max(counts)

        # Store High Count (converted to dr)
        high_count *= ((rad[-1] - rad[0]) / num_rad)

        if t_i == 0:
            r_widths0[i] = high_count
        elif t_i == 1:
            r_widths1[i] = high_count
        elif t_i == 2:
            r_widths2[i] = high_count
        elif t_i == 3:
            r_widths5[i] = high_count
        else:
            r_widths10[i] = high_count

    ####### Radii #######
    rs[i] = peak_rad



## Use These Frames ##
rate = 5 # 5 works better, but is very slow
start = 10
max_frame = util.find_max_frame()
frame_range = np.array(range(start, max_frame + 1, rate))

# Store These Values
thresholds = [0.05, 0.1, 0.2, 0.5, 1.0]

az_widths0 = mp_array("d", len(frame_range)) # 0.05
az_widths1 = mp_array("d", len(frame_range)) # 0.1
az_widths2 = mp_array("d", len(frame_range)) # 0.2
az_widths5 = mp_array("d", len(frame_range)) # 0.5
az_widths10 = mp_array("d", len(frame_range)) # 1.0

r_widths0 = mp_array("d", len(frame_range)) # 0.05
r_widths1 = mp_array("d", len(frame_range)) # 0.1
r_widths2 = mp_array("d", len(frame_range)) # 0.2
r_widths5 = mp_array("d", len(frame_range)) # 0.5
r_widths10 = mp_array("d", len(frame_range)) # 1.0

rs = mp_array("d", len(frame_range)) # Track 'r' for aspect ratio

# Call Function Over Time
pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

p = Pool(10)
p.map(get_extents, pool_args)
p.terminate()

##### PLOTTING #####

# Plot Parameters
linewidth = 4
fontsize = 14

my_dpi = 100
alpha = 0.5

def make_azimuthal_extent_plot():
    # Figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    # Curves
    plot.plot(frame_range, az_widths0, linewidth = linewidth, label = "0.05")
    plot.plot(frame_range, az_widths1, linewidth = linewidth, label = "0.1")
    plot.plot(frame_range, az_widths2, linewidth = linewidth, label = "0.2")
    plot.plot(frame_range, az_widths5, linewidth = linewidth, label = "0.5")
    plot.plot(frame_range, az_widths10, linewidth = linewidth, label = "1.0")

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Azimuthal Extents (degrees)", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    plot.legend(loc = "upper left")

    # Limits
    plot.xlim(0, frame_range[-1])
    #plot.ylim(0.0, 1.0)

    # Save + Close
    plot.savefig("azimuthal_extents.png")
    plot.show()

    plot.close()

def make_radial_extent_plot():
    # Figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    # Curves
    plot.plot(frame_range, r_widths0, linewidth = linewidth, label = "0.05")
    plot.plot(frame_range, r_widths1, linewidth = linewidth, label = "0.1")
    plot.plot(frame_range, r_widths2, linewidth = linewidth, label = "0.2")
    plot.plot(frame_range, r_widths5, linewidth = linewidth, label = "0.5")
    plot.plot(frame_range, r_widths10, linewidth = linewidth, label = "1.0")

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Radial Extents", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    plot.legend(loc = "upper left")

    # Limits
    plot.xlim(0, frame_range[-1])
    #plot.ylim(0.0, 1.0)

    # Save + Close
    plot.savefig("radial_extents.png")
    plot.show()

    plot.close()

def make_aspect_ratio_plot():
    # Figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    # Collect Data
    aspect_ratios0 = np.nan_to_num(np.array(rs) * np.array(az_widths0) / np.array(r_widths0))
    aspect_ratios1 = np.nan_to_num(np.array(rs) * np.array(az_widths1) / np.array(r_widths1))
    aspect_ratios2 = np.nan_to_num(np.array(rs) * np.array(az_widths2) / np.array(r_widths2))
    aspect_ratios5 = np.nan_to_num(np.array(rs) * np.array(az_widths5) / np.array(r_widths5))
    aspect_ratios10 = np.nan_to_num(np.array(rs) * np.array(az_widths10) / np.array(r_widths10))

    # Curves
    plot.plot(frame_range, aspect_ratios0, linewidth = linewidth, label = "0.05")
    plot.plot(frame_range, aspect_ratios1, linewidth = linewidth, label = "0.1")
    plot.plot(frame_range, aspect_ratios2, linewidth = linewidth, label = "0.2")
    plot.plot(frame_range, aspect_ratios5, linewidth = linewidth, label = "0.5")
    plot.plot(frame_range, aspect_ratios10, linewidth = linewidth, label = "1.0")

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Aspect Ratios", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    plot.legend(loc = "upper left")

    # Limits
    plot.xlim(0, frame_range[-1])
    #plot.ylim(0.0, 1.0)

    # Save + Close
    plot.savefig("aspect_ratios.png")
    plot.show()

    plot.close()


make_azimuthal_extent_plot()
make_radial_extent_plot()
make_aspect_ratio_plot()

