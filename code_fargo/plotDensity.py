"""
plots density over time
then, makes movies
"""

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import fromfile

import subprocess

# Plot Parameters
my_dpi = 100
cmap = "RdYlBu_r"
clim = [10000000, -10000000] # Set later

# Get FARGO Parameters
num_frames = 480
num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

# Find minimum and maximum
for i in range(num_frames):
    density = np.log(fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
    minimum = np.min(density)
    maximum = np.max(density)
    if minimum < clim[0]:
        clim[0] = minimum
    if maximum > clim[1]:
        clim[1] = maximum
print clim

# Plot Each File
for i in range(num_frames):
    density = np.log(fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))

    ### (1b.i) Magnitude in Polar ###
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    ax = fig.add_subplot(111, polar = True)
    subset_r = range(0, int(num_rad / 1.3))
    subset_t = range(0, int(num_theta))
    result = ax.pcolormesh(theta[subset_t], rad[subset_r], (density[subset_r])[:, subset_t], cmap = cmap)
    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    fontsize = 14

    plot.title("Gas Density at Timestep %d" % i, fontsize = fontsize + 1)
    plot.savefig("gasDensity/frame%03d_polar.png" % i, bbox_inches = 'tight', dpi = my_dpi)
    #plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)

    ### (1b.ii) Magnitude in "r-theta Cartesian" ###
    fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
    ax = fig.add_subplot(111)
    subset_r = range(0, int(num_rad / 1.3))
    subset_t = range(0, int(num_theta))
    result = ax.pcolormesh(theta[subset_t], rad[subset_r], (density[subset_r])[:, subset_t], cmap = cmap)
    fig.colorbar(result)
    result.set_clim(clim[0], clim[1])

    fontsize = 14

    plot.title("Gas Density at Timestep %d" % i, fontsize = fontsize + 1)
    plot.savefig("gasDensity/frame%03d_cart.png" % i, bbox_inches = 'tight', dpi = my_dpi)
    #plot.show()
    plot.close(fig) # Close Figure (to avoid too many figures)


#### Make Movies ####
command = "avconv -framerate 40 -f image2 -vf scale=-2:720 -i gasDensity/frame%03d_polar.png -b 65536k gasDensity/gasDensity_polar.mov"
split_command = command.split()
subprocess.Popen(split_command)

command = "avconv -framerate 40 -f image2 -vf scale=-2:720 -i gasDensity/frame%03d_cart.png -b 65536k gasDensity/gasDensity_cart.mov"
split_command = command.split()
subprocess.Popen(split_command)
