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

def find_peak(averagedProfile):
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    peak_rad_outer_index = np.argmax(averagedProfile[outer_disk_start:])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_value = averagedProfile[peak_index]

    return peak_rad, peak_index, peak_value

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

def find_azimuthal_peak(azimuthalProfile):
    ### Use Smoothed Profile ###
    twenty_degrees = (np.pi / 180.0) * 20.0
    kernel_size = np.searchsorted(rad, twenty_degrees) # kernel corresponds to 20 degrees

    smoothed_profile = smooth(azimuthalProfile, kernel_size)
    peak_theta_index = np.argmax(smoothed_profile)

    peak_theta = theta[peak_theta_index]
    peak_value = azimuthalProfile[peak_theta_index]

    return peak_theta, peak_theta_index, peak_value

#### Data ####

def get_data(frame):
    # Find Peak in Radial Profile (in Outer Disk)
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    peak_rad, peak_index, peak_density = find_peak(averagedDensity)
    min_rad, min_density = find_min(averagedDensity, peak_rad)

    peak_theta, peak_theta_index, peak_azimuthal_density = find_azimuthal_peak(density[peak_index])

    if len(sys.argv) > 2:
        # Supply central theta as an argument
        peak_theta = float(sys.argv[2])

    # Gather Azimuthal Profiles
    pressure = density * (scale_height)**2 * np.power(rad, -1)[:, None]

    num_profiles = 5
    spread = 30.0 # half-width

    radial_theta = np.linspace(peak_theta - spread, peak_theta + spread, num_profiles)
    radial_theta[radial_theta < 0] += 360.0
    radial_theta[radial_theta > 360] -= 360.0

    radial_indices = [np.searchsorted(theta, this_theta * (np.pi / 180.0)) for this_theta in radial_theta]
    radial_profiles = [pressure[:, radial_index] for radial_index in radial_indices]

    return radial_theta, radial_profiles

##### PLOTTING #####

# Make Directory
directory = "radialPressure"
try:
    os.mkdir(directory)
except:
    print "Directory Already Exists"

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

alpha = 0.65
fontsize = 14
linewidth = 4

def make_plot(frame, radial_theta, radial_profiles, show = False):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    ### Plot ###
    x_min = 1.0
    x_max = 2.5
    x_min_i = np.searchsorted(rad, x_min)
    x_max_i = np.searchsorted(rad, x_max)

    max_pressure = 0

    for this_theta, radial_profile in zip(radial_theta, radial_profiles):
        print this_theta
        plot.plot(rad, radial_profile, linewidth = linewidth, alpha = alpha, label = "%.1f" % this_theta)

        # Store max for ylim
        this_max_pressure = np.max(radial_profile[x_min_i : x_max_i])
        if this_max_pressure > max_pressure:
            max_pressure = this_max_pressure


    # Axis
    angles = np.linspace(0, 2 * np.pi, 7)
    degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

    plot.xlim(x_min, x_max)
    plot.ylim(0, 1.1 * max_pressure)

    # Annotate
    this_title = readTitle()
    plot.xlabel("Radius", fontsize = fontsize + 2)
    plot.ylabel("Pressure", fontsize = fontsize)
    plot.title("Orbit %d: %s" % (orbit, this_title), fontsize = fontsize + 1)

    plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 1.0)) # outside of plot)

    # Save and Close
    plot.savefig("%s/radial_pressure_%04d.png" % (directory, frame), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    if frame_number == -1:
        # Plot Sample
        max_frame = util.find_max_frame()
        sample = np.linspace(10, max_frame, 10) # 10 evenly spaced frames
        for i in sample:
            radial_theta, radial_profiles = get_data(i)
            make_plot(i, radial_theta, radial_profiles)
    else:
        # Plot Single
        radial_theta, radial_profiles = get_data(frame_number)
        make_plot(frame_number, radial_theta, radial_profiles, show = True)
else:
    # Search for maximum frame
    density_files = glob.glob("gasdens*.dat")
    max_frame = find_max_frame()
    num_frames = max_frame + 1

    #for i in range(num_frames):
    #    make_plot(i)

    #### ADD TRY + CATCH BLOCK HERE!!!!! ####

    #p = Pool() # default number of processes is multiprocessing.cpu_count()
    #p.map(make_plot, range(num_frames))
    #p.terminate()

    #### Make Movies ####
    #make_movies()
