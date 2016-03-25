"""
Usage:
python plotSampleAveragedDensity.py
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

### Choose Sample ###
sample_name = "by200"
sample_range = np.linspace(0, 2000, 11)

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

##### PLOTTING #####

# Make Directory
directory = "averagedDensity"
try:
    os.mkdir(directory)
except:
    print "Directory Already Exists"

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 4

def add_to_plot(frame, show = False):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Axis
    x = rad
    prefix = ""
    plot.xlim(float(fargo_par["Rmin"]) - 0.05, 1.0)
    plot.ylim(0, 5.0)
        
    # Data
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    log_density = np.log(fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    averagedDensity = np.average(density, axis = 1) / surface_density

    ### Plot ###
    plot.plot(x, averagedDensity, linewidth = linewidth, label = "%d" % frame)

    # Annotate
    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel("Azimuthally Averaged Density", fontsize = fontsize)
    #plot.title("", fontsize = fontsize + 1)

def make_plot(prefix = "", show = False):
    # Annotate
    plot.legend()

    # Save and Close
    plot.savefig("%s/%savg_density_%s.png" % (directory, prefix, sample_name), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()


##### Plot One File or All Files #####

# Set up figure


### Plot Sample ###
# Full Range
fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
for i in sample_range:
    add_to_plot(i)

make_plot(prefix = "inner_", show = True)