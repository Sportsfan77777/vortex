"""
adds up negative vorticity (ideally in the vortex) over time

Usage:
python plotTotalVorticityOverTime.py
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool
from multiprocessing import Array as mp_array

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
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

def sum_vorticity(args):
    i, frame = args

    # Get Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero
    vrad = np.array(fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
    vtheta = np.array(fromfile("gasvtheta%d.dat" % frame).reshape(num_rad, num_theta))

    # Get Background Data
    #vtheta_keplerian = np.array(fromfile("gasvtheta0.dat").reshape(num_rad, num_theta))

    # Subtract off Keplerian velocity (and rotate back into non-rotating frame???)
    #vtheta -= (vtheta_keplerian)
    vorticity = util.velocity_curl(vrad, vtheta, rad, theta, frame = 0)

    # Add up vorticity
    dr = rad[1] - rad[0] # assumes arithmetic grid
    d_phi = theta[1] - theta[0]

    radial_vorticity = np.average(vorticity, axis = 1)
    total_vorticity = np.sum((dr * d_phi) * vortex_rad[:, None] * radial_vorticity)
    
    # Print Update
    print "%d: %.4f" % (frame, total_vorticity)

    # Store Data
    vorticity_over_time[i] = total_vorticity


## Use These Frames ##
rate = 10 # 5 works better, but is very slow
start = 10
max_frame = util.find_max_frame()
frame_range = np.array(range(start, max_frame + 1, rate))

# Gather Here (with multiprocessing)
vorticity_over_time = mp_array("d", len(frame_range))

pool_args = [(i, frame) for i, frame in enumerate(frame_range)]

#for i, frame in enumerate(frame_range):
#   sum_vorticity((i, frame))

p = Pool(10)
p.map(sum_vorticity, pool_args)
p.terminate()

## Pickle to combine later ##

pickle.dump(np.array(frame_range), open("integrated_disk_vorticity_frames_taper%d.p" % taper_time, "wb"))
pickle.dump(np.array(vorticity_over_time), open("integrated_disk_vorticity_values_taper%d.p" % taper_time, "wb"))

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
    plot.plot(frame_range, vorticity_over_time, linewidth = linewidth)

    # Annotate
    this_title = readTitle()
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Integrated Disk Vorticity", fontsize = fontsize)
    plot.title(this_title, fontsize = fontsize)

    #plot.legend(loc = "upper right")

    # Limits
    plot.xlim(0, frame_range[-1])
    #plot.ylim(0.0, 1.0)

    # Save + Close
    plot.savefig("integratedDiskVorticity.png")
    plot.show()

    plot.close()


make_plot()



