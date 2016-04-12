"""
Compares tapering "accretion rate" to usual measurements (involving viscosity or radial velocity)

Usage:
python diagnoseAccretion.py
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
from readTitle import readTitle

## Set file names ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    planet_fn = "planet1.dat"
else:
    # fargo
    planet_fn = "planet0.dat"

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

mass_taper = int(round(float((fargo_par["MassTaper"])), 0))

##### PLOTTING #####

# Plot Parameters
my_dpi = 100

fontsize = 14
linewidth = 4

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)

    ### Data ###
    xs = range(mass_taper)

    # Measure 1: Diff Mass Accretion
    data = np.loadtxt(planet_fn)
    masses = (data[:, 5]) # Planet Mass from 0 to Mass Taper
    tapering_accretion = np.diff(masses[:mass_taper + 1])

    print masses[mass_taper + 1], masses[-1]

    # Measure 2: Viscosity
    viscous_accretion = 3.0 * np.pi * surface_density * viscosity

    # Measure 3: Radial Velocity
    near_planet = 1.05
    radius_near_planet = np.searchsorted(rad, near_planet)
    
    rate = 10
    xs_vrad = range(0, mass_taper, rate)
    radial_velocity_at_planet = np.zeros(len(xs_vrad))

    for i, frame in enumerate(xs_vrad):
        vrad = (fromfile("gasvrad%d.dat" % frame).reshape(num_rad, num_theta))
        radial_velocity_at_planet[i] = np.average(vrad[radius_near_planet, :])

    radius = 1.0
    radial_accretion = 2.0 * np.pi * radius * surface_density * radial_velocity_at_planet

    # Curves
    plot.plot(xs, tapering_accretion, color = "blue", label = "taper", linewidth = linewidth)
    plot.plot([xs[0], xs[-1]], [viscous_accretion, viscous_accretion], color = "black", label = "viscous", linewidth = linewidth)
    plot.plot(xs_vrad, radial_velocity_at_planet, color = "red", label = "v_rad (+)", linewidth = linewidth) # Positive
    plot.plot(xs_vrad, -radial_velocity_at_planet, color = "orange", label = "v_rad (-)", linewidth = linewidth) # Negative

    # Limits
    plot.xlim(xs[0], xs[-1])
    plot.yscale('log')

    # Annotate
    this_title = readTitle()
    plot.xlabel("Time", fontsize = fontsize)
    plot.ylabel("Accretion Rate", fontsize = fontsize)
    plot.title("Accretion for %s" % (this_title), fontsize = fontsize + 1)

    plot.legend(loc = "upper right", bbox_to_anchor = (1.36, 1.0)) # outside of plot

    # Save + Close
    plot.savefig("accretion_diagnosis.png", bbox_inches = 'tight', dpi = my_dpi)
    plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

make_plot()