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

def get_concentration_times(pressure_gradients, width = 4 * scale_height, omega = 1, dust_density = 1.2, stokes_number = 1.0):
    stokes_factor = 1.0 / (stokes_number + stokes_number**(-1))
    drift_velocities = stokes_factor * pressure_gradients / dust_density / omega

    return width / drift_velocities

#### Data ####

def get_data(frame):
    # Read Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density
    averagedDensity = np.average(density, axis = 1)

    # Limit to Outer Disk
    rad_planet = np.searchsorted(rad, 1.0)
    outer_rad = rad[rad_planet:]
    outer_density = density[rad_planet:]

    # Find Peak
    avg_outer_density = np.average(outer_density, axis = 1)
    peak_index = np.argmax(avg_outer_density)
    peak_rad = np.searchsorted(outer_rad, peak_index)

    # Only Look for ArgMax-es near peak
    scale_height_i = np.searchsorted(rad, rad[0] + scale_height)

    mask = np.ones(np.shape(outer_density))
    mask[peak_index - scale_height_i : peak_index + scale_height_i] = 0 # only look within +/- 1 scale height

    # Mask Pressure and Get ArgMax-es
    pressure = outer_density * (scale_height)**2 * np.power(outer_rad, -1)[:, None]
    masked_pressure = np.ma.array(pressure, mask = mask)

    maxima = np.ma.argmax(masked_pressure, axis = 0, fill_value = 0)

    maxima_i = maxima
    theta_i = np.indices(np.shape(maxima))

    # Take Pressure Gradients
    num_h = 4
    half_width_i = num_h * scale_height_i
    half_width = num_h * scale_height
    inner_pressure_gradients = np.abs(pressure[maxima_i, theta_i] - pressure[maxima_i - half_width_i, theta_i]) / half_width
    outer_pressure_gradients = np.abs(pressure[maxima_i, theta_i] - pressure[maxima_i + half_width_i, theta_i]) / half_width

    inner_pressure_gradients = inner_pressure_gradients[0]
    outer_pressure_gradients = outer_pressure_gradients[0]

    # Convert to Concentration Time
    omega = peak_rad**(-1.5)
    inner_concentration_times = get_concentration_times(nner_pressure_gradients, width = half_width, omega = omega)
    outer_concentration_times = get_concentration_times(outer_pressure_gradients, width = half_width, omega = omega)

    return inner_concentration_times, outer_concentration_times

##### PLOTTING #####

# Make Directory
directory = "dustConcentrations"
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

def make_plot(frame, inner_pressure_gradients, outer_pressure_gradients, show = False):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    ### Plot ###
    x = np.linspace(0, 360, num_theta)
    plot.plot(x, inner_pressure_gradients, linewidth = linewidth, label = "Inner")
    plot.plot(x, outer_pressure_gradients, linewidth = linewidth, label = "Outer")

    # Axis
    angles = np.linspace(0, 360, 7)
    degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

    plot.xticks(angles)
    plot.xlim(0, 360)
    plot.ylim(0, 1.1 * np.max(inner_pressure_gradients))

    # Annotate
    this_title = readTitle()
    plot.xlabel(r"$\phi$", fontsize = fontsize + 2)
    plot.ylabel("Dust Concentration Time (in planet orbits)", fontsize = fontsize)
    plot.title("Orbit %d: %s" % (orbit, this_title), fontsize = fontsize + 1)

    plot.legend(loc = "upper right") # outside of plot)

    # Save and Close
    plot.savefig("%s/dust_concentration_%04d.png" % (directory, frame), bbox_inches = 'tight', dpi = my_dpi)
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
            inner_pressure_gradients, outer_pressure_gradients = get_data(i)
            make_plot(i, inner_pressure_gradients, outer_pressure_gradients)
    else:
        # Plot Single
        inner_pressure_gradients, outer_pressure_gradients = get_data(frame_number)
        make_plot(frame_number, inner_pressure_gradients, outer_pressure_gradients, show = True)
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

