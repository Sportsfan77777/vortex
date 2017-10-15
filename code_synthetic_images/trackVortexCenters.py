"""
Usage: trackVortexCenters.py

Plots vortex centers over time.
"""

import pickle
import glob
from multiprocessing import Pool

import numpy as np

import argparse
from pylab import fromfile

import util

directories = ["cm", "hcm", "mm", "hmm", "hum", "um"]
sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001])

### Input Parameters ###

def new_argument_parser(description = "Generate input for synthetic images."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    return parser

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get FARGO Parameters ###
directory = "../%s-size" % directories[0]
fargo_par = util.get_pickled_parameters(directory = directory) # Retrieve parameters from *.par file

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

scale_height = fargo_par["AspectRatio"]

### Get Input Parameters ###

# Frames
if len(args.frames) == 1:
    frame_range = args.frames
elif len(args.frames) == 3:
    start = args.frames[0]; end = args.frames[1]; rate = args.frames[2]
    frame_range = range(start, end + 1, rate)

# Number of Cores 
num_cores = args.num_cores

### Analyis Output ##

centers_m1 = np.zeros((len(frame_range), len(directories)))
centers_m2 = np.zeros((len(frame_range), len(directories)))

### Task Functions ###

def retrieve_density(frame, size):
    """ Step 0: Retrieve density """

    fn_i = "../%s-size/gasddens%d.dat" % (size, frame)
    density = fromfile(fn_i).reshape(num_rad, num_theta)

    return density

def method1(density):
    """ argmax """
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    argmax = np.argmax(density_segment)
    arg_r, arg_phi = np.unravel_index(argmax, np.shape(density_segment))

    return theta[arg_phi]

def method2(density):
    """ center of threshold """
    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    # Get peak in azimuthal profile
    avg_density = np.average(density_segment, axis = 1) # avg over theta
    arg_peak = np.argmax(avg_density)
    peak_rad = np.searchsorted(rad, 1.1 + arg_peak)

    # Average over half a scale height
    half_width = 0.25 * scale_height
    zoom_start = np.searchsorted(peak_rad - half_width)
    zoom_end = np.searchsorted(peak_rad + half_width)

    # Zoom in on peak
    density_sliver = density[zoom_start : zoom_end]
    avg_density_sliver = np.average(density_sliver, axis = 0) # avg over rad

    # Center vortex (loosely)
    middle = np.searchsorted(theta, np.pi)
    arg_center = np.argmax(avg_density_sliver)

    shift_i = int(middle - center_i)
    avg_density_sliver = np.roll(avg_density_sliver, shift_i, axis = 1)

    # Pick threshold
    threshold = np.percentile(density_sliver, 95)

    left_i = np.searchsorted(avg_density_sliver, side = "left")
    right_i = np.searchsorted(avg_density_sliver, side = "right")

    middle_i = (left_i + right_i) / 2.0

    return theta[middle_i] + theta[shift_i]

def full_procedure(index):
    directory = "../%s-size" % directories[index]

    for i, frame in enumerate(frame_range):
        density = retrieve_density(frame, directories[index])

        center_m1 = method1(density)
        center_m2 = method2(density)

        centers_m1[i, index] = center_m1
        centers_m2[i, index] = center_m2


###############################################################################

##### PLOTTING #####

fontsize = 14
linewidth = 3

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

def make_plot(centers_mi, method, show = True):

    # Set up figure
    fig = plot.figure()

    ### Line Plot ###
    for i, directory in enumerate(directories):
        plot.plot(frame_range, centers_mi[:, i], c = colors[i], linewidth = linewidth, label = "directory")

    # Annotate Axes
    plot.xlabel("Time (planet orbits)", fontsize = fontsize)
    plot.ylabel(r"$\phi$", fontsize = fontsize)
    if method == 3:
        plot.title("Difference in Vortex Centers", fontsize = fontsize + 1)
    else:
        plot.title("Vortex Centers: Method %d" % method, fontsize = fontsize + 1)

    plot.legend(loc = "upper right")

    # Axes
    plot.xlim(frame_range[0], frame_range[-1])

    if method == 3:
        plot.ylim(-180, 180)

        angles = np.linspace(-180, 180, 7)
        plot.yticks(angles)
    else:
        plot.ylim(0, 360)

        angles = np.linspace(0, 360, 7)
        plot.yticks(angles)

    # Save, Show, and Close
    save_fn = "vortexCenters_method%d.png" % method
    plot.savefig(save_fn, bbox_inches = 'tight')

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


### Generate Synthetic Input ###

# Iterate through frames

if num_cores > 1:
    indices = range(6)

    p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
    p.map(full_procedure, indices)
    p.terminate()
else:
    for i, directory in enumerate(directories):
        full_procedure(i)

make_plot(centers_m1, 1)
make_plot(centers_m2, 2)
make_plot(centers_m2 - centers_m1, 3)
