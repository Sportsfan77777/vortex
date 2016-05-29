"""
plots excess mass in difference images over time

Usage:
python plotExcessMassOverTime.py
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

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

planet_mass = float(fargo_par["PlanetMass"])
mass_taper = float(fargo_par["MassTaper"])

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

def get_excess_mass(args):
    # Unwrap Args
    i, frame = args

    # Get Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero
    background_density = (fromfile("gasdens%d.dat" % (frame - 1)).reshape(num_rad, num_theta)) / surface_density_zero

    diff_density = density - background_density
    diff_density[diff_density < 0] = 0 # only include excess

    # Extract Near Vortex
    averagedDensity = np.average(density, axis = 1)
    peak_rad, peak_density = find_peak(averagedDensity)

    vortex_start = np.max([1.0, peak_rad - 5.0 * scale_height])
    vortex_end = peak_rad + 5.0 * scale_height

    vortex_start_i = np.searchsorted(rad, vortex_start)
    vortex_end_i = np.searchsorted(rad, vortex_end)

    vortex_rad = rad[vortex_start_i : vortex_end_i]
    vortex_diff_density = diff_density[vortex_start_i : vortex_end_i]

    vortex_excess = np.average(vortex_diff_density, axis = 1)

    # Add up mass
    dr = rad[1] - rad[0] # assumes arithmetic grid
    d_phi = theta[1] - theta[0]

    excess_mass = np.sum((dr * d_phi) * vortex_rad[:, None] * vortex_diff_density)

    # Get Peak
    peak_diff_density = np.max(vortex_excess)

    # Print Update
    print "%d: %.4f, %.4f" % (frame, excess_mass, peak_diff_density)

    # Store Data
    mass_over_time[i] = excess_mass
    peak_over_time[i] = peak_diff_density

def record_triggers():
    triggers = {}

    thresholds = [0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.20]
    for threshold in thresholds:
        frame_i = np.searchsorted(mass_over_time, threshold)
        triggers[threshold] = frame_range[frame_i]

    # Pickle
    pickle.dump(triggers, open("triggers.p", "wb"))

## Use These Frames ##
profile_function = lambda x : (planet_mass / 0.001) * np.power(np.sin((np.pi / 2) * mass_taper), 2)
profile = [profile_function(x) for x in range(mass_taper)]

rate = 1
start = np.searchsorted(profile, 0.08) # start with 0.08 Jupiter masses
end = np.searchsorted(profile, 0.52) # end with 0.52 Jupiter masses
frame_range = np.array(range(start, end, rate))

#mass_over_time = np.zeros(len(frame_range))
#peak_over_time = np.zeros(len(frame_range))

mass_over_time = mp_array("d", len(frame_range))
peak_over_time = mp_array("d", len(frame_range))

#for i, frame in enumerate(frame_range):
#    get_excess_mass((i, frame))

pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

p = Pool(10)
p.map(get_excess_mass, pool_args)
p.terminate()

max_mass = np.max(mass_over_time)
max_peak = np.max(peak_over_time)

## Mark Triggers ##

record_triggers()

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
    plot.plot(frame_range, mass_over_time, linewidth = linewidth)
    #plot.plot(frame_range, peak_over_time, linewidth = linewidth - 1, label = "Peak")

    # Reference Lines
    #plot.plot([0, frame_range[-1]], 0.10 * max_mass * np.ones(2), linewidth = 2, color = "black")
    #plot.plot([0, frame_range[-1]], 0.10 * max_peak * np.ones(2), linewidth = 1, color = "black")

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Excess Mass", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    #plot.legend(loc = "upper right")

    # Limits
    plot.xlim(0, frame_range[-1])
    #plot.ylim(0.0, 1.0)

    # Save + Close
    plot.savefig("triggerExcessMass.png")
    plot.show()

    plot.close()


make_plot()


