"""
plots azimuthal max as a function of radius near the vortex

python trackSpiralWave.py
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
from readTitle import readTitle

save_directory = "sprialWaves"

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

### Helper Functions ###
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter

def getWaveLocation(density, radii, start_radius = 1.10, end_radius = 2.30):
    # Find Limits
    start_i = np.searchsorted(radii, start_radius)
    end_i = np.searchsorted(radii, end_radius)

    radii = radii[start_i : end_i] # truncate to selected range

    # Find maxes
    density_near_vortex = density[start_i : end_i]
    wave_locations = np.zeros(len(density_near_vortex))

    guess_theta = 2 * np.pi
    delta_theta = 10.0 * (np.pi / 180.0)
    for i, r_i in enumerate(radii):
        # Guess Bounds
        upper_guess = guess_theta # + delta_theta
        lower_guess = guess_theta - delta_theta # search only below previous location

        # Mask Values Away From Wave 
        density_ring = density_near_vortex[i]
        mask = np.ones(len(density_ring)) # 1 is masked

        if lower_guess < 0.0:
            lower_guess_i = np.searchsorted(theta, lower_guess + 2 * np.pi)
            upper_guess_i = np.searchsorted(theta, upper_guess)

            mask[lower_guess_i : ] = 0
            mask[ : upper_guess_i] = 0

        elif upper_guess > 2 * np.pi:
            lower_guess_i = np.searchsorted(theta, lower_guess)
            upper_guess_i = np.searchsorted(theta, upper_guess - 2 * np.pi)

            mask[lower_guess_i : ] = 0
            mask[ : upper_guess_i] = 0

        else:
            lower_guess_i = np.searchsorted(theta, lower_guess)
            upper_guess_i = np.searchsorted(theta, upper_guess)
            
            mask[lower_guess_i : upper_guess_i] = 0

        masked_density_ring = np.ma.array(density_ring, mask = mask)
        wave_index = np.ma.argmax(masked_density_ring, fill_value = -1.0)

        # Store Result + Update Next Guess
        wave_locations[i] = theta[wave_index]
        guess_theta = wave_locations[i]

    return radii, wave_locations

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
linewidth = 4
fontsize = 14
alpha = 0.5

my_dpi = 100

def make_plot(frame, radii, wave_locations, show = False):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    # Data
    #kernel_size = 5
    #smoothed_wave_locations = smooth(wave_locations, kernel_size)

    # Curves
    plot.plot(radii, wave_locations, linewidth = linewidth)

    # Axis
    plot.xlim(radii[0], radii[-1])

    angles = np.linspace(0, 2 * np.pi, 7)
    degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

    plot.ylim(0, 2 * np.pi)
    plot.yticks(angles, degree_angles)

    # Annotate
    this_title = readTitle()
    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)
    plot.title("Sprial Wave Location at Orbit %d\n%s" % (orbit, this_title), fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("%s/sprialWave_%04d.png" % (save_directory, frame), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    if frame_number == -1:
        # Plot Sample
        max_frame = util.find_max_frame()
        sample = np.linspace(50, max_frame, 10) # 10 evenly spaced frames
        for i in sample:
            density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta)) / surface_density_zero
            radii, wave_locations = util.getWaveLocation(density, rad)
            make_plot(i, density, radii, wave_locations)
    else:
        # Plot Single
        density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta)) / surface_density_zero
        radii, wave_locations = util.getWaveLocation(density, rad)
        make_plot(frame_number, density, radii, wave_locations, show = True)
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



