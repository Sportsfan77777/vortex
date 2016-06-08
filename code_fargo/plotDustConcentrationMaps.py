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
from scipy.optimize import curve_fit

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import rcParams # replace with rc ???
from pylab import fromfile

import util
from readTitle import readTitle

import colormaps as cmaps
plot.register_cmap(name = 'viridis', cmap = cmaps.viridis)
plot.register_cmap(name = 'inferno', cmap = cmaps.inferno)
plot.register_cmap(name = 'plasma', cmap = cmaps.plasma)
plot.register_cmap(name = 'magma', cmap = cmaps.magma)

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

def calculate_stokes_number(size, dust_density, gas_surface_density):
    return (np.pi / 2) * (size * dust_density) / gas_surface_density

def get_concentration_time(pressure_gradient, gas_surface_density, omega, width = scale_height, stokes_number = 1.0):
    stokes_factor = 1.0 / (stokes_number + stokes_number**(-1))
    drift_velocities = stokes_factor * pressure_gradient / gas_surface_density / omega[:, None]

    return width / drift_velocities

#### Data ####

def get_data(frame):
    # Read Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    # Convert to pressure and pressure gradient
    dr = rad[1] - rad[0]

    sound_speeds = (scale_height)**2 * np.power(rad, -1)
    pressure = density * sound_speeds[:, None]
    pressure_gradient = np.abs(np.diff(pressure, axis = 0)) / dr

    # Get Concentration Times
    omega = np.power(rad, -1.5)
    concentration_times = get_concentration_time(pressure_gradient, density[1:], omega[1:])

    return concentration_times

## Array ##

##### PLOTTING #####

# Make Directory
directory = "dustConcentrationMaps"
try:
    os.mkdir(directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "inferno"
clim = [1, 4]

rcParams['figure.figsize'] = 5, 10
my_dpi = 100

alpha = 0.65
fontsize = 14
linewidth = 4

def make_plot(frame, concentration_times, show = False):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    ### Plot ###
    result = plot.pcolormesh(rad[1:], theta, np.transpose(np.log(concentration_times)), cmap = cmap)

    cbar = fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    # Axis
    plot.xlim(1.0, 2.5)
    
    angles = np.linspace(0, 2 * np.pi, 7)
    degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

    plot.ylim(0, 2 * np.pi)
    plot.yticks(angles, degree_angles)

    # Annotate
    this_title = readTitle()
    plot.xlabel(r"$\phi$", fontsize = fontsize + 2)
    plot.ylabel("Dust Concentration Maps", fontsize = fontsize)
    plot.title("Orbit %d: %s" % (orbit, this_title), fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("%s/dust_concentration_map_%04d.png" % (directory, frame), bbox_inches = 'tight', dpi = my_dpi)
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
            concentration_times = get_data(i)
            make_plot(i, concentration_times)
    else:
        # Plot Single
        concentration_times = get_data(frame_number)
        make_plot(frame_number, concentration_times, show = True)
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

