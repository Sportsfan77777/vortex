"""
plot 2-D density maps

python plotBothDensityMaps.py
python plotBothDensityMaps.py frame_number
python plotBothDensityMaps.py -1 <<<===== Plots a sample
python plotBothDensityMaps.py -m
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np


import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
from readTitle import readTitle

save_directory = "bothMeasuredViscosityMaps"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1, 0]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

# Keplerian Velocity Profiles
vtheta0 = (fromfile("gasvtheta0.dat").reshape(num_rad, num_theta))
dust_vtheta0 = (fromfile("gasdvtheta0.dat").reshape(num_rad, num_theta))

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "RdYlBu_r"
gas_clim = [-9, -4]
dust_clim = [-9, -4]

fontsize = 14
my_dpi = 100


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
        angles = np.linspace(0, 2 * np.pi, 7)
        degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

        plot.ylim(0, 2 * np.pi)
        plot.yticks(angles, degree_angles)
        if axis == "zoom":
            x = (rad - 1) / scale_height
            prefix = "zoom_"
            plot.xlim(0, 40) # to match the ApJL paper
            xlabel = r"($r - r_p$) $/$ $h$"
        else:
            x = rad
            prefix = ""
            plot.xlim(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]))
            xlabel = "Radius"

        # Data
        density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
        normalized_density = density / surface_density_zero

        vrad = (fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
        vtheta = (fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta)) - vtheta0

        # Calculate Viscosity Parameter ()
        measured_viscosity = vrad * vtheta
        log_measured_viscosity = np.log10(measured_viscosity)

        ### Plot ###
        result = plot.pcolormesh(x, theta, np.transpose(log_measured_viscosity), cmap = cmap)
        fig.colorbar(result)
        result.set_clim(gas_clim[0], gas_clim[1])

        # Annotate
        this_title = readTitle()
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel(r"$\phi$", fontsize = fontsize)
        plot.title("Measured Viscosity Map at Orbit %d\n%s" % (orbit, this_title), fontsize = fontsize + 1)

        ##### Dust Density #####
        plot.subplot(1, 2, 2)

        # Axis
        angles = np.linspace(0, 2 * np.pi, 7)
        degree_angles = ["%d" % d_a for d_a in np.linspace(0, 360, 7)]

        plot.ylim(0, 2 * np.pi)
        plot.yticks(angles, degree_angles)
        if axis == "zoom":
            x = (rad - 1) / scale_height
            prefix = "zoom_"
            plot.xlim(0, 40) # to match the ApJL paper
            xlabel = r"($r - r_p$) $/$ $h$"
        else:
            x = rad
            prefix = ""
            plot.xlim(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]))
            xlabel = "Radius"

        # Data
        dust_density = (fromfile("gasddens%d.dat" % i).reshape(num_rad, num_theta))
        normalized_dust_density = dust_density / surface_density_zero

        dust_vrad = (fromfile("gasdvrad%d.dat" % i).reshape(num_rad, num_theta))
        dust_vtheta = (fromfile("gasdvtheta%d.dat" % i).reshape(num_rad, num_theta)) - dust_vtheta0

        # Calculate Viscosity Parameter ()
        measured_dust_viscosity = dust_vrad * dust_vtheta
        log_measured_dust_viscosity = np.log10(measured_dust_viscosity)

        ### Plot ###
        result = plot.pcolormesh(x, theta, np.transpose(log_measured_dust_viscosity), cmap = cmap)
        fig.colorbar(result)
        result.set_clim(dust_clim[0], dust_clim[1])

        # Annotate
        this_title = readTitle()
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel(r"$\phi$", fontsize = fontsize)
        plot.title("Measured Dust Viscosity Map at Orbit %d\n%s" % (orbit, this_title), fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/%sbothMeasuredViscosityMaps_%04d.png" % (save_directory, prefix, i), bbox_inches = 'tight', dpi = my_dpi)
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
        sample = np.linspace(50, max_frame, 10) # 10 evenly spaced frames
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



