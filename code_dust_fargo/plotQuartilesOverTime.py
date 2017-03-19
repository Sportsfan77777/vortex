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

def get_quartiles(args):
    # Unwrap Args
    i, frame = args

    # Get Dust Data
    density = (fromfile("gasddens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero
    background_density = (fromfile("gasddens0.dat").reshape(num_rad, num_theta)) / surface_density_zero

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
    vortex_background_density = background_density[vortex_start_i : vortex_end_i]

    # Look for 2x initial density (mask everything else)
    not_overdense = np.array(vortex_density < 2.0 * vortex_background_density)
    masked_density = (np.ma.array(vortex_density, mask = not_overdense)).compressed() # return non-masked data only

    # Get Data (if not all masked)
    if not np.all(not_overdense):
        m50 = np.percentile(masked_density, 50) # Note: mask has no effect on np.percentile
        m75 = np.percentile(masked_density, 75)
        m90 = np.percentile(masked_density, 90)
        m99 = np.percentile(masked_density, 99)

        # Store Data
        track50[i] = m50
        track75[i] = m75
        track90[i] = m90
        track99[i] = m99

        # Print Update
        print "%d: %.4f, %.4f, %.4f, %.4f" % (frame, m50, m75, m90, m99)

## Use These Frames ##
rate = 5 # 5 works better, but is very slow
start = 10
max_frame = util.find_max_frame()
frame_range = np.array(range(start, max_frame + 1, rate))

track50 = mp_array("d", len(frame_range))
track75 = mp_array("d", len(frame_range))
track90 = mp_array("d", len(frame_range))
track99 = mp_array("d", len(frame_range))

pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

p = Pool(10)
p.map(get_quartiles, pool_args)
p.terminate()

##### PLOTTING #####

# Plot Parameters
linewidth = 4
fontsize = 14

my_dpi = 100
alpha = 0.5

def make_plot():
    # Figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    # Curves
    plot.plot(frame_range, track99, linewidth = linewidth, label = "99")
    plot.plot(frame_range, track90, linewidth = linewidth, label = "90")
    plot.plot(frame_range, track75, linewidth = linewidth, label = "75")
    plot.plot(frame_range, track50, linewidth = linewidth, label = "50")

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Dust Density", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    plot.legend(loc = "upper left")

    # Limits
    plot.xlim(0, frame_range[-1])
    #plot.ylim(0.0, 1.0)

    # Save + Close
    plot.savefig("quartiles.png")
    plot.show()

    plot.close()


make_plot()


