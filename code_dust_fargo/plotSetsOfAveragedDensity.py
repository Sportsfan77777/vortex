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

save_directory = "averagedDensity"

name = "averagedDensities%d.p"
ids = ["5524", "6420"]

### Load Pickle Files ###

# Master Entry #
master = pickle.load(open(name % ids[0], "rb"))
frame_range = master["frame_range"]
rad = master["rad"]
scale_height = master["scale_height"]

# All #
cubes = []
for i, id_i in enumerate(ids):
    entry = pickle.load(open(name % ids[i], "rb"))
    cube_i = entry["cube"]
    cubes.append(cube_i)

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 4

def make_plot(i, frame, show = False):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, frame, axis):
        # Orbit Number
        #time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
        #orbit = int(round(time / (2 * np.pi), 0)) * i
        orbit = frame

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
            plot.xlim(0, x[-1])
            xlabel = "Radius"

        ### Plot ###
            
        # Data
        for cube in cubes:
            averagedDensity = cube[i]
            plot.plot(x, averagedDensity, linewidth = linewidth)        

        # Annotate
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel("Azimuthally-Averaged\n Normalized Density", fontsize = fontsize)
        plot.title(r"$t = %d$ $\rm{orbits}" % (orbit), fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/%savg_density_%04d.png" % (save_directory, prefix, frame), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    # Two Plots
    choose_axis(i, frame, "normal")
    choose_axis(i, frame, "zoom")

##### Plot One File Per Frame #####

for i, frame in enumerate(frame_range):
    make_plot(i, frame)


