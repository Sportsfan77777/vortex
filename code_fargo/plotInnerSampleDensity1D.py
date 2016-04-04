"""
For FARGO2D1D, plot 1D-density profile.

Usage:
python plotSampleDensity1D.py
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

from readTitle import readTitle

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

rad = np.loadtxt("used_rad1D.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

##### PLOTTING #####

# Make Directory
directory = "density1D"
try:
    os.mkdir(directory)
except:
    print "Directory Already Exists"

# Plot Parameters
rcParams['figure.figsize'] = 5, 10
my_dpi = 100

fontsize = 14
linewidth = 4

def add_to_plot(frame, choice = "normal", show = False):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Axis
    x = rad
    prefix = ""
    plot.xlim(float(fargo_par["Rmin1D"]) - 0.05, 1.0)
    plot.ylim(0, 5.0)
        
    # Data
    averagedDensity = (np.loadtxt("gasdens.ascii_rad.%d.dat" % frame))[:, 1] / surface_density

    ### Plot ###
    plot.plot(x, averagedDensity, linewidth = linewidth, label = "%d" % frame)

    # Annotate
    this_title = readTitle()
    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel("1-D Density", fontsize = fontsize)
    plot.title("%s" % this_title, fontsize = fontsize + 1)

def make_plot(prefix = "", show = False):
    # Annotate
    plot.legend()

    # Save and Close
    plot.savefig("%s/%sdensity1D_%s.png" % (directory, prefix, sample_name), bbox_inches = 'tight', dpi = my_dpi)
    if show:
        plot.show()


##### Plot One File or All Files #####

### Plot Sample ###
# Full Range
fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
for i in sample_range:
    add_to_plot(i)

make_plot(prefix = "inner_", show = True)