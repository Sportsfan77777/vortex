"""
plot 2-D density AND excess vorticity maps

python plotVorticity.py
python plotVorticity.py frame_number
python plotVorticity.py -1  <<<===== Plots a sample
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

import util
from readTitle import readTitle

save_directory = "diffComboMaps"

## Check frame ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    ref_frame = 0
else:
    # fargo
    ref_frame = 1

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

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
d_cmap = "afmhot"
d_clim = [0, 2]

v_cmap = "afmhot_r"
v_clim = [-0.2, 0.0]

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
        #gs = gridspec.GridSpec(1, 2)
        #ax1 = fig.add_subplot(gs[0])
        #ax2 = fig.add_subplot(gs[1])

        ###### Density ######
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
        density = (fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta)) / surface_density_zero
        background_density = (fromfile("gasdens%d.dat" % (i - 1)).reshape(num_rad, num_theta)) / surface_density_zero

        diff_density = density - background_density

        ### Plot ###
        result = plot.pcolormesh(x, theta, np.transpose(diff_density), cmap = d_cmap)
        fig.colorbar(result)
        result.set_clim(d_clim[0], d_clim[1])

        # Annotate
        this_title = readTitle()
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel(r"$\phi$", fontsize = fontsize)
        plot.title("Gas Density Difference Map at Orbit %d\n%s" % (orbit, this_title), fontsize = fontsize + 1)

        ###### Excess Vortensity ######
        # Select Plot
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
        vtheta_keplerian = np.array(fromfile("gasvtheta0.dat").reshape(num_rad, num_theta))

        ### FRAME 1 ###
        
        vrad = np.array(fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta))
        vtheta = np.array(fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta))

        # Subtract off Keplerian velocity (and rotate back into non-rotating frame???)
        vtheta -= (vtheta_keplerian)
        vorticity = util.velocity_curl(vrad, vtheta, rad, theta, frame = 0)

        # Divide out Angular Frequency (for Rossby Number)
        vorticity = vorticity / (np.array([r**(-1.5) for r in rad[:-1]]))[:, None]

        ### FRAME 2 ###

        background_vrad = np.array(fromfile("gasvrad%d.dat" % (i - 1)).reshape(num_rad, num_theta))
        background_vtheta = np.array(fromfile("gasvtheta%d.dat" % (i - 1)).reshape(num_rad, num_theta))

        # Subtract off Keplerian velocity (and rotate back into non-rotating frame???)
        background_vtheta -= (vtheta_keplerian)
        background_vorticity = util.velocity_curl(background_vrad, background_vtheta, rad, theta, frame = 0)

        # Divide out Angular Frequency (for Rossby Number)
        background_vorticity = background_vorticity / (np.array([r**(-1.5) for r in rad[:-1]]))[:, None]

        # Remove (vorticity > 0) from background
        background_vorticity[background_vorticity > 0] = 0

        ### TAKE DIFFERENCE
        diff_vorticity = vorticity - background_vorticity

        ### Plot ###
        result = plot.pcolormesh(x, theta, np.transpose(diff_vorticity), cmap = v_cmap)
    
        fig.colorbar(result)
        result.set_clim(v_clim[0], v_clim[1])

        # Annotate
        plot.xlabel(xlabel, fontsize = fontsize)
        plot.ylabel(r"$\phi$", fontsize = fontsize)
        plot.title("Vorticity (v - vK) Map at Orbit %d\n%s" % (orbit, this_title), fontsize = fontsize + 1)

        # Save and Close
        plot.savefig("%s/%sdiffComboMap_%03d.png" % (save_directory, prefix, i), bbox_inches = 'tight', dpi = my_dpi)
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
        sample = np.linspace(1, max_frame, 125) # 125 evenly spaced frames
        
        #for i in sample:
        #    make_plot(i)

        p = Pool(5) # default number of processes is multiprocessing.cpu_count()
        p.map(make_plot, sample)
        p.terminate()
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


