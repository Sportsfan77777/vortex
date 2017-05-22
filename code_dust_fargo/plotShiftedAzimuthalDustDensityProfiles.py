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

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

rad = np.linspace(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]), float(fargo_par["Nrad"]))
num_rad = len(rad)

num_theta = float(fargo_par["Nsec"])
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

def get_data(frame, size):
    # Find Peak in Radial Profile (in Outer Disk)
    density = (fromfile("shifted_gasddens%d_%s.npy" % (frame, size))[:-10].reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    peak_rad, peak_density = find_peak(averagedDensity)
    min_rad, min_density = find_min(averagedDensity, peak_rad)

    if (len(sys.argv) > 3):
        peak_rad = float(sys.argv[3])

    # Gather Azimuthal Profiles
    num_profiles = 7
    spread = 1.0 * scale_height # half-width

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = [density[azimuthal_index, :] for azimuthal_index in azimuthal_indices]

    return azimuthal_radii, azimuthal_profiles

##### PLOTTING #####

# Make Directory
directory = "azimuthalDustDensity"
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

def make_plot(frame, azimuthal_radii, azimuthal_profiles, show = False):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    ### Plot ###
    for radius, azimuthal_profile in zip(azimuthal_radii, azimuthal_profiles):
        plot.plot(theta, azimuthal_profile, linewidth = linewidth, alpha = alpha, label = "%.3f" % radius)

    # Axis
    angles = np.linspace(0, 2 * np.pi, 7)
    degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

    plot.xlim(0, 2 * np.pi)
    plot.xticks(angles, degree_angles)

    # Annotate
    this_title = "Size: %s" % size #readTitle()
    plot.xlabel(r"$\phi$", fontsize = fontsize + 2)
    plot.ylabel("Shifted Azimuthal Density", fontsize = fontsize)
    plot.title("Orbit %d: %s" % (orbit, this_title), fontsize = fontsize + 1)

    plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 1.0)) # outside of plot)

    # Save and Close
    plot.savefig("%s/shifted_azimuthal_density_%04d_%s.png" % (directory, frame, size), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    size = sys.argv[2]
    if frame_number == -1:
        # Plot Sample
        max_frame = util.find_max_frame()
        sample = np.linspace(10, max_frame, 10) # 10 evenly spaced frames
        for i in sample:
            azimuthal_radii, azimuthal_profiles = get_data(i, size)
            make_plot(i, azimuthal_radii, azimuthal_profiles)
    else:
        # Plot Single
        azimuthal_radii, azimuthal_profiles = get_data(frame_number, size)
        make_plot(frame_number, azimuthal_radii, azimuthal_profiles, show = True)
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
