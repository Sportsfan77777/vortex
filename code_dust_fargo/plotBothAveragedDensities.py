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

save_directory = "bothAveragedDensities"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1, 0]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

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

def make_plot(frame, show = False):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, axis):
        # Orbit Number
        time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
        orbit = int(round(time / (2 * np.pi), 0)) * i

        # Set up figure
        fig = plot.figure(figsize = (1400 / my_dpi, 600 / my_dpi), dpi = my_dpi)
        
        ### Gas Density ###
        # Select Plot
        plot.subplot(1, 2, 1)

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
            
        # Data
        density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
        averagedDensity = np.average(density, axis = 1) / surface_density

        ### Plot ###
        plot.plot(x, averagedDensity, linewidth = linewidth)

        # Annotate
        this_title = readTitle()
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel("Azimuthally-Averaged\n Gas Density", fontsize = fontsize)
        plot.title(r"$t = %d$ $\rm{orbits}\n%s" % (orbit, this_title), fontsize = fontsize + 1)

        ### Dust Density ###
        # Select Plot
        plot.subplot(1, 2, 2)
        
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
            
        # Data
        density = (fromfile("gasddens%d.dat" % i).reshape(num_rad, num_theta))
        averagedDensity = np.average(density, axis = 1) / surface_density

        ### Plot ###
        plot.plot(x, averagedDensity, linewidth = linewidth)

        # Annotate
        this_title = readTitle()
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel("Azimuthally-Averaged\n Dust Density", fontsize = fontsize)
        plot.title(r"$t = %d$ $\rm{orbits}\n%s" % (orbit, this_title), fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/%sboth_avg_density_%04d.png" % (save_directory, prefix, i), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    choose_axis(i, "normal")
    choose_axis(i, "zoom")

##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    if frame_number == -1:
        # Plot Sample
        max_frame = util.find_max_frame()
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
