"""
The mass accretion rate should be related to the radial velocity of the flow
2pi sigma r v_r = 3pi sigma viscosity

(2/3) r v_r ?= viscosity
Let's find out! 

Usage:
python diagnoseViscosity.py frame
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
viscosity = float(fargo_par["Viscosity"])

##### PLOTTING #####

# Make Directory
directory = "diagnosingViscosity"
try:
    os.mkdir(directory)
except:
    print "Directory Already Exists"

# Plot Parameters
my_dpi = 100

fontsize = 14
linewidth = 4

def make_plot(frame):
    # Orbit Number
    time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
    orbit = int(round(time / (2 * np.pi), 0)) * frame

    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    # Data
    fargo_fn = "fargo2D1D"
    if os.path.exists(fargo_fn):
        # fargo2D1D
        vrad = fromfile("gasvrad1D%d.dat" % frame)
        pseudo_viscosity = (2.0 / 3) * (np.multiply(rad, vrad)) # should be equal to viscosity
        avg_pseudo_viscosity = np.average(pseudo_viscosity, axis = 1) # radial pseudo-viscosity
    else:
        # fargo
        vrad = (fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
        pseudo_viscosity = (2.0 / 3) * (np.multiply(rad[:, None], vrad)) # should be equal to viscosity
        avg_pseudo_viscosity = np.average(pseudo_viscosity, axis = 1) # radial pseudo-viscosity

    # Curves
    plot.plot(rad, avg_pseudo_viscosity, color = "blue", label = "+", linewidth = linewidth) # positive
    plot.plot(rad, -avg_pseudo_viscosity, color = "red", label = "-", linewidth = linewidth) # negative
    plot.plot([rad[0], rad[-1]], [viscosity, viscosity], color = "black", linewidth = linewidth - 1)

    # Limits
    plot.xlim(rad[0], rad[-1])
    plot.yscale('log')

    # Annotate
    plot.xlabel("Radius", fontsize = fontsize)
    plot.ylabel("Pseudo Viscosity", fontsize = fontsize)
    plot.title("Orbit %d" % orbit, fontsize = fontsize + 1)

    plot.legend(loc = "upper right")

    # Save + Close
    plot.savefig("%s/viscosity_diagnosis_%04d.png" % (directory, frame), bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


frame = int(sys.argv[1])
make_plot(frame)