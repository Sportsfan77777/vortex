"""
plot azimuthally averaged density
then, makes movies
"""

import os
import subprocess

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import rcParams # replace with rc ???
from pylab import fromfile

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 4

# Get FARGO Parameters
num_frames = 390 # Calculate this instead using glob search?

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density = 6.366 * 10**(-4)
scale_height = 0.06 # Find this in .par file?

# Make Directory
directory = "averagedDensity"
try:
    os.mkdir(directory)
except:
    print "Directory Already Exists"

def make_plot(frame):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, axis):
        # Orbit Number
        orbit = 2 * i # This could be different for different parameters!!!!!!!

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
            
        # Data
        density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
        log_density = np.log(fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
        averagedDensity = np.average(density, axis = 1) / surface_density

        ### Plot ###
        plot.plot(x, averagedDensity, linewidth = linewidth)

        # Annotate
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel("Azimuthally Averaged Density", fontsize = fontsize)
        plot.title("Averaged Gas Density at Orbit %d" % orbit, fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("averagedDensity/%savg_density_%03d.png" % (prefix, i), bbox_inches = 'tight', dpi = my_dpi)
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    choose_axis(i, "normal")
    choose_axis(i, "zoom")

# Plot Each File
for i in range(num_frames):
    make_plot(i)

#### Make Movies ####
# Movie Parameters
fps = 40

path = "averagedDensity/avg_density_%03d.png"
output = "averagedDensity/averagedDensity.mov"

path = "averagedDensity/zoom_avg_density_%03d.png"
zoom_output = "averagedDensity/averagedDensity.mov"

# Movie Command
command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
split_command = command.split()
subprocess.Popen(split_command)

command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, zoom_path, zoom_output)
split_command = command.split()
subprocess.Popen(split_command)

