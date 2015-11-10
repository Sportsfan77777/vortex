"""
plot azimuthally averaged density
then, makes movies

Usage:
python plotAveragedDensity.py frame_number <== plot one frame
python plotAveragedDensity.py -m <== make a movie instead of plots (plots already exist)
python plotAveragedDensity.py <== plot all frames and make a movie

"""

import sys
import os
import subprocess
import glob

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import rcParams # replace with rc ???
from pylab import fromfile

### Movie Commands ###
def make_movies():
    # Movie Parameters
    fps = 40

    path = "averagedDensity/avg_density_%03d.png"
    output = "averagedDensity/averagedDensity.mov"

    zoom_path = "averagedDensity/zoom_avg_density_%03d.png"
    zoom_output = "averagedDensity/averagedDensity_zoom.mov"

    # Movie Command
    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    split_command = command.split()
    subprocess.Popen(split_command)

    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, zoom_path, zoom_output)
    split_command = command.split()
    subprocess.Popen(split_command)

# Make only movies and then return
if (len(sys.argv) > 2) and (sys.argv[2] == "-m"):
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
fargo_par = pickle.load(open(params_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

# Search for maximum frame
density_files = glob.glob("gasdens*.dat")
max_frame = 0
for d_f in density_files:
    name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
    frame_number = int(d_f[7:]) # just 999
    if frame_number > max_frame:
        max_frame = frame_number

num_frames = max_frame # Calculate this instead using glob search?

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

def make_plot(frame):
    # For each frame, make two plots (one with normal 'r' and one with '(r - 1) / h')
    def choose_axis(i, axis):
        # Orbit Number
        orbit = int(round(fargo_par["Ninterm"] * fargo_par["DT"] / (2 * np.pi), 0))

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

##### Plot One File or All Files #####

if len(sys.argv) > 2:
    frame_number = sys.argv[2]
    make_plot(frame_number)
else:
    for i in range(num_frames):
        make_plot(i)

    #### Make Movies ####
    make_movies()

