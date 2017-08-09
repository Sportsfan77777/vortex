"""
plot 2-D density maps

python plotDensityMaps.py
python plotDensityMaps.py frame_number
python plotDensityMaps.py -1 <<<===== Plots a sample
python plotDensityMaps.py -m
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

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

save_directory = "shiftedDustDensityMaps"

# Units
radius = 5.0
radius_unit = radius * (1.496 * 10**13) # (AU / cm)

# Vortex Locations
# T_growth = 1000
centers = {}
centers["cm"] = 119.5
centers["hcm"] = 214.9
centers["mm"] = 180.2
centers["hmm"] = 184.5
centers["hum"] = 185.6
centers["um"] = 184.5

# T_growth = 10
centers["cm"] = 120.9
centers["hcm"] = 80.2
centers["mm"] = 55.3
centers["hmm"] = 55.9
centers["hum"] = 58.2
centers["um"] = 56.7

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

rad = np.loadtxt("radial.dat") / radius_unit
theta = np.loadtxt("azimuthal.dat")

num_rad = len(rad)
num_theta = len(theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

### Helper Functions ###
middle = np.searchsorted(theta, np.pi)

def find_argmax(density):
    # Returns Azimuthal Argmax
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    argmax = np.argmax(density_segment)
    arg_r, arg_phi = np.unravel_index(argmax, np.shape(density_segment))

    return arg_phi

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "RdYlBu_r"
clim = [0, 2]

fontsize = 14
my_dpi = 100


def make_plot(frame, show = False):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, axis):
        # Orbit Number
        time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
        orbit = int(round(time / (2 * np.pi), 0)) * i

        # Set up figure
        fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
        ax = fig.add_subplot(111)

        # Axis
        angles = np.linspace(0, 2 * np.pi, 7)
        degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

        plot.ylim(0, 2 * np.pi)
        plot.yticks(angles, degree_angles)
        if axis == "zoom":
            x = (rad - 1) / scale_height
            prefix = "zoom_"
            plot.xlim(0, 40) # to match the ApJL paper
            xlabel = r"($r - r_p$) $/$ $h$"
        else:
            x = rad
            prefix = ""
            plot.xlim(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]))
            xlabel = "Radius"

        # Data
        density = (fromfile("gasdens%d_%s.dat" % (i, size)).reshape(num_rad, num_theta))
        normalized_density = density / surface_density_zero

        # Rotate Vortex to Center
        #location_i = find_argmax(normalized_density)
        location_phi = (np.pi / 180.0) * centers[size]
        location_i = np.searchsorted(theta, location_phi)

        shift_i = int(middle - location_i)
        normalized_density = np.roll(normalized_density, shift_i, axis = 1)

        ### Plot ###
        result = ax.pcolormesh(x, theta, np.transpose(normalized_density), cmap = cmap)
        fig.colorbar(result)
        result.set_clim(clim[0], clim[1])

        # Annotate
        this_title = readTitle()
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel(r"$\phi$", fontsize = fontsize)
        plot.title("Shifted %s-size Gas Density Map at Orbit %d\n%s" % (size, orbit, this_title), fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/%sshiftedGasDensityMap%04d_%s.png" % (save_directory, prefix, i, size), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

        # Save Shifted Data
        #### Write this part! ####
        shifted_data = np.roll(density, shift_i, axis = 1)
        shift_savename = "shifted_gasdens%d_%s.p" % (i, size)
        #np.save(shift_savename, shifted_data)
        pickle.dump(shifted_data, open(shift_savename, 'r'))

    i = frame
    choose_axis(i, "normal")
    choose_axis(i, "zoom")


##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    size = sys.argv[2]
    if frame_number == -1:
        # Plot Sample
        max_frame = util.find_max_frame()
        sample = np.linspace(50, max_frame, 10) # 10 evenly spaced frames
        for i in sample:
            make_plot(i)
    else:
        # Plot Single
        make_plot(frame_number, show = True)
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



