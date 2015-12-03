"""
plot 2-D vorticity maps (really vortensity = vorticity / density)

python plotVorticity.py
python plotVorticity.py frame_number
python plotVorticity.py -m
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


save_directory = "vorticityMaps"

### Movie Commands ###
def make_movies():
    # Movie Parameters
    fps = 5

    path = save_directory + "/vorticityMap_%03d.png"
    output = save_directory + "/vorticityMap.mov"

    zoom_path = save_directory + "/zoom_vorticityMap_%03d.png"
    zoom_output = save_directory + "/vorticityMap_zoom.mov"

    # Movie Command
    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    split_command = command.split()
    subprocess.Popen(split_command)

    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, zoom_path, zoom_output)
    split_command = command.split()
    subprocess.Popen(split_command)

# Make only movies and then return
if (len(sys.argv) > 1) and (sys.argv[1] == "-m"):
    make_movies()
    # Terminate
    quit()


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

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

# Curl function
def curl(v_rad, v_theta, rad, theta):
    """ z-component of the curl (because this is a 2-D simulation)"""
    ### Start Differentials ###
    d_rad = np.diff(rad)
    d_theta = np.diff(theta)

    dv_rad = np.diff(v_rad, axis = 1)
    dv_theta = np.diff(rad[:, None] * v_theta, axis = 0)
    ### End Differentials ###

    # z-Determinant
    partial_one = dv_theta / d_rad[:, None]
    partial_two = dv_rad / d_theta

    z_curl = (partial_one[:, 1:] - partial_two[1:, :]) / rad[1:, None]
    return z_curl


# Curl function
def old_curl(v_rad, v_theta, rad, theta):
    """ z-component of the curl (because this is a 2-D simulation)"""
    ### Start Differentials ###
    # d_r
    d_rad = np.diff(rad)

    # d_t
    d_theta = np.diff(theta)

    order = 1

    # dv_rad
    dv_rad = np.diff(v_rad, axis = 1, n = order)

    # dv_theta
    dv_theta = np.diff(rad * v_theta, axis = 0, n = order)

    ### End Differentials ###

    # z-Determinant
    partial_one = dv_theta / d_rad[:, None] # Note: dr is one shorter than rad
    partial_two = dv_rad / d_theta # Note: dt is d_theta, not d_time!!!!!!!!!
    #partial_one = rad[:-1, None] * dv_theta[:-1] / dr # Note: dr is one shorter than rad
    #partial_two = dv_rad[:-1, :] / dt # Note: dt is d_theta, not d_time!!!!!!!!!

    print "One"
    print np.median(partial_one), partial_one
    print "Two"
    print np.median(partial_two), partial_two

    # Source: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
    z_curl = (partial_one[:, 1:] - partial_two[1:, :]) / rad[1:]

    print "Curl"
    print np.median(z_curl), z_curl

    return z_curl


##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "RdYlBu_r"
clim = [-2, 2]

fontsize = 14
my_dpi = 100

def make_plot(frame):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, axis):
        # Orbit Number
        time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
        orbit = int(round(time / (2 * np.pi), 0)) * i

        # Set up figure
        fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
        ax = fig.add_subplot(111)

        # Axis
        if axis == "zoom":
            x = (rad - 1) / scale_height
            prefix = "zoom_"
            plot.xlim(0, 40) # to match the ApJL paper
            plot.ylim(0, 2 * np.pi)
            xlabel = r"($r - r_p$) $/$ $h$"
        else:
            x = rad
            prefix = ""
            plot.xlim(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]))
            plot.ylim(0, 2 * np.pi)
            xlabel = "Radius"

        # Data
        density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
        normalized_density = density / surface_density_zero

        vrad = (fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
        vtheta = (fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta))

        w = curl(vrad, vtheta, rad, theta)

        print len(w[0,:]), len(w[:,0])
        print len(normalized_density[0,:]), len(normalized_density[:,0])

        print normalized_density[:,0]

        ### Plot ###
        #result = ax.pcolormesh(x, theta, np.transpose(w), cmap = cmap)
        result = ax.pcolormesh(x, theta, np.transpose(w / normalized_density[:len(w[0,:]), :len(w[:,0])]), cmap = cmap)
        #result = ax.pcolormesh(x, theta, np.transpose(w / normalized_density), cmap = cmap)
        fig.colorbar(result)
        result.set_clim(clim[0], clim[1])

        # Annotate
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel(r"$\phi$", fontsize = fontsize)
        plot.title("Gas Density Map at Orbit %d" % orbit, fontsize = fontsize + 1)

        # Save and Close
        #plot.savefig("%s/%svorticityMap_%03d.png" % (save_directory, prefix, i), bbox_inches = 'tight', dpi = my_dpi)
        plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    #choose_axis(i, "normal")
    choose_axis(i, "zoom")



##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    make_plot(frame_number)
else:
    # Search for maximum frame
    density_files = glob.glob("gasdens*.dat")
    max_frame = 0
    for d_f in density_files:
        name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
        frame_number = int(name[7:]) # just 999
        if frame_number > max_frame:
            max_frame = frame_number
    num_frames = max_frame + 1

    #for i in range(num_frames):
    #    make_plot(i)

    p = Pool() # default number of processes is multiprocessing.cpu_count()
    p.map(make_plot, range(num_frames))
    p.terminate()

    #### Make Movies ####
    make_movies()


