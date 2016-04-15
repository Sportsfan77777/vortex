"""
plot azimuthally averaged density
then, makes movies

Usage:
python plotAveragedDensity.py frame_number <== plot one frame
python plotAveragedDensity.py -m <== make a movie instead of plots (plots already exist)
python plotAveragedDensity.py -1 <<<===== Plots a sample
python plotAveragedDensity.py <== plot all frames and make a movie

"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool
import time

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import rcParams # replace with rc ???
from pylab import fromfile

from readTitle import readTitle

name = "RWI_Trigger"

directories = []
frame_numbers = []
save_directory = "./"

# Store Current Directory for Later
working_directory = os.getcwd()

# Keep track of xs and ys, and titles as labels
xs = []
ys = []
titles = []

for (directory, frame) in zip(directories, frame_numbers):
    os.chdir(directory)

    ### Get FARGO Parameters ###
    # Create param file if it doesn't already exist
    param_fn = "params.p"
    if not os.path.exists(param_fn):
        command = "python pickleParameters.py"
        split_command = command.split()
        subprocess.Popen(split_command)
        time.delay(2) # delay two seconds to give time for command to execute (to pickle parameters)
    fargo_par = pickle.load(open(param_fn, "rb"))

    num_rad = np.loadtxt("dims.dat")[-2]
    num_theta = np.loadtxt("dims.dat")[-1]

    rad = np.loadtxt("used_rad.dat")[:-1]
    theta = np.linspace(0, 2 * np.pi, num_theta)

    surface_density = float(fargo_par["Sigma0"])
    scale_height = float(fargo_par["AspectRatio"])

    # Get Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    averagedDensity = np.average(density, axis = 1) / surface_density

    this_title = readTitle()

    # Add Data to Collective
    xs.append(rad)
    ys.append(averaged_density)
    titles.append(this_title)

##### PLOTTING #####

os.chdir(working_directory)

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 4

def make_plot(show = False):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, axis):
        # Orbit Number
        time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
        orbit = int(round(time / (2 * np.pi), 0)) * i

        # Set up figure
        fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

        # Axis
        if axis == "zoom":
            x = (rad - 1) / scale_height
            prefix = "zoom_"
            plot.xlim(0, 40) # to match the ApJL paper
            plot.ylim(0, 1.2)
            xlabel = r"($r - r_p$) $/$ $h$"
        else:
            x = rad
            prefix = ""
            xlabel = "Radius"

        ### Plot ###
        for (this_x, this_y, frame, this_title) in zip(xs, ys, frame_numbers, titles):
            label = "%s for %d" % (this_title, frame)
            plot.plot(this_x, this_y, label = this_title, linewidth = linewidth)

        # Annotate
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel("Azimuthally Averaged Density", fontsize = fontsize)
        plot.title("%s" % (name), fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%savg_density_%s.png" % (prefix, name), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    choose_axis(i, "normal")
    choose_axis(i, "zoom")

##### Plot One File or All Files #####

make_plot(frame_number, show = True)
