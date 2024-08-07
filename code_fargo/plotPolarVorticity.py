"""
plot 2-D vorticity maps (really vortensity = vorticity / density)

python plotPolarVorticity.py
python plotPolarVorticity.py frame_number
python plotPolarVorticity.py -1  <<<===== Plots a sample
python plotPolarVorticity.py -m
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

## Check frame ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    frame = 0
else:
    # fargo
    frame = 1

save_directory = "polarVortensityMaps"

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

res = [int(fargo_par["Nrad"]), int(fargo_par["Nsec"])]
boundary = fargo_par["InnerBoundary"]

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "RdYlBu_r"
clim = [0, 0.2]

fontsize = 14
my_dpi = 100

def make_plot(frame, show = False):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, axis):
        # Orbit Number
        time = float(fargo_par["Ninterm"]) * float(fargo_par["DT"])
        orbit = int(round(time / (2 * np.pi), 0)) * i

        # Set up figure
        fig = plot.figure(figsize = (700 / my_dpi, 600 / my_dpi), dpi = my_dpi)
        ax = fig.add_subplot(111, polar = True)

        # Axis
        if axis == "zoom":
            prefix = "zoom_"
            rmax = 2.4 # to match ApJL paper
        else:
            prefix = ""
            rmax = float(fargo_par["Rmax"])

        plot.xticks([], []) # Angles
        plot.yticks([rmax], ["%.1f" % rmax]) # Max Radius
        plot.ylim(0, rmax)

        # Data
        density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta))
        normalized_density = density / surface_density_zero

        vrad = (fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
        vtheta = (fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta))

        vorticity = util.velocity_curl(vrad, vtheta, rad, theta) / normalized_density[1:, 1:]

        ### Plot ###
        result = ax.pcolormesh(theta, rad, vorticity, cmap = cmap)
    
        fig.colorbar(result)
        result.set_clim(clim[0], clim[1])

        # Annotate
        this_title = readTitle()
        plot.title("Vortensity Map at Orbit %d\n%s" % (orbit, this_title), fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/%svortensityMap_%03d.png" % (save_directory, prefix, i), bbox_inches = 'tight', dpi = my_dpi)
        if show:
            plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

    i = frame
    choose_axis(i, "normal")
    choose_axis(i, "zoom")



##### Plot One File or All Files #####

def find_max_frame():
    max_frame = 0
    for d_f in density_files:
        name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
        frame_number = int(name[7:]) # just 999
        if frame_number > max_frame:
            max_frame = frame_number
    return max_frame

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    if frame_number == -1:
        # Plot Sample
        max_frame = find_max_frame()
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

    p = Pool() # default number of processes is multiprocessing.cpu_count()
    p.map(make_plot, range(num_frames))
    p.terminate()

    #### Make Movies ####
    #make_movies()


