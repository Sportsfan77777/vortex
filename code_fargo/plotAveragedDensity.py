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
num_frames = 480
num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

# Make Directory
directory = "averagedDensity"
try:
    os.mkdir(snapshot_dir)
except:
    print "Directory Already Exists"

# Plot Each File
for i in range(num_frames):
    density = np.log(fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
    averagedDensity = np.average(density, axis = 1)

    ### Plot ###
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    plot.plot(rad, averagedDensity, linewidth = linewidth)

    # Annotate
    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel("Azimuthally Averaged Density", fontsize = fontsize)
    plot.title("Averaged Gas Density at Timestep %d" % i, fontsize = fontsize + 1)

    # Save and Close
    plot.savefig("averagedDensity/avg_density_%03d.png" % i, bbox_inches = 'tight', dpi = my_dpi)
    #plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


#### Make Movies ####
command = "avconv -framerate 40 -f image2 -vf scale=-2:720 -i averagedDensity/avg_density_%03d.png -b 65536k gasDensity/gasDensity_polar.mov"
split_command = command.split()
subprocess.Popen(split_command)

