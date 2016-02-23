"""
plots contributions to torque for each cell

Usage:
python plotTorqueMap.py
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

save_directory = "torqueMaps"

# Constants #
BigG = 1
planet_mass = 0.005

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


# Torque Function Element by Element
def torque(radius, theta, density):
    """ 'torque' = 'r cross F' """
    # Differentials
    d_rad = np.diff(rad)
    d_theta = np.diff(theta)

    # Relevant Vectors + Quantities
    r_element = np.array([np.outer(radius, np.cos(theta)), np.outer(radius, np.sin(theta))]) # to fluid element 
    r_diff = np.array([np.outer(radius, np.cos(theta)) - 1, np.outer(radius, np.sin(theta))]) # to fluid element 

    dist_sq = np.linalg.norm(r_diff)**2

    print np.shape(r_element), np.shape(r_diff)

    print np.shape(dist_sq)

    # Torque
    coeff = BigG * planet_mass * density / dist_sq
    direction = np.cross(r_element, r_diff, axis = 0)

    torque_density = coeff * direction
    #area = d_rad**2 * d_theta

    return torque_density

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 3

cmap = "RdYlBu_r"
clim = [-16, 0]
#clim = [-2, 2] # direction-only clim

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
        density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
        normalized_density = density / surface_density_zero

        log_torque = np.log(np.abs(torque(rad, theta, normalized_density)))
        print np.shape(log_torque)

        ### Plot ###
        result = ax.pcolormesh(x, theta, np.transpose(log_torque), cmap = cmap)
    
        fig.colorbar(result)
        result.set_clim(clim[0], clim[1])

        # Annotate
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel(r"$\phi$", fontsize = fontsize)
        plot.title("Torque Map at Orbit %d" % orbit, fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/%storqueMap_%03d.png" % (save_directory, prefix, i), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    choose_axis(i, "normal")
    choose_axis(i, "zoom")

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    if frame_number == -1:
        # Plot Sample
        max_frame = find_max_frame()
        sample = np.linspace(10, max_frame, 10) # 10 evenly spaced frames
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

