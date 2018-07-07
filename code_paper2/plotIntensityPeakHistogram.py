"""
plot 2-D density maps

python plotDensityMaps.py
python plotDensityMaps.py frame_number
python plotDensityMaps.py -1 <<<===== Plots a sample
python plotDensityMaps.py -m
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import gridspec

from pylab import rcParams
from pylab import fromfile

import util
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "intensityPeakHistogram",
                         help = 'save directory (default: intensityPeakHistogram)')

    # Save Data
    parser.add_argument('--data', dest = "save_data", action = 'store_true', default = False,
                         help = 'save data or not (default: do not save)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'max_y for plot (default: None)')

    parser.add_argument('--retrieve', dest = "measure", action = 'store_false', default = True,
                         help = 'measure or retrieve peak offset (default: measure)')
    parser.add_argument('-t', dest = "threshold", type = float, default = 0.5,
                         help = 'threshold for marking edges of vortex with intensity (default: 0.5)')
    parser.add_argument('--print', dest = "print_histogram", action = 'store_true', default = False,
                         help = 'measure or retrieve peak offset (default: measure)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'labelsize of plot annotations (default: 16)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'linewidths in plot (default: 3)')
    parser.add_argument('--alpha', dest = "alpha", type = float, default = 0.65,
                         help = 'line transparency in plot (default: 0.65)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get ID%04d Parameters ###
fn = "id%04d_par.p" % (args.id_number)
fargo_par = pickle.load(open(fn, "rb"))

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
dust_surface_density_zero = surface_density_zero / 100
disk_mass = 2 * np.pi * dust_surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

planet_radius = fargo_par["Radius"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Save Data
save_data = args.save_data

# Plot Parameters (variable)
show = args.show
max_y = args.max_y
#if max_y is None:
#    max_y = 1

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

measure = args.measure
threshold = args.threshold
print_histogram = args.print_histogram

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Data ###

def get_peak_offset(args):
    """return peak in azimuthal profile of vortex"""
    # Unpack
    i, frame = args

    intensity_polar = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number)
    arg_r, arg_phi = az.get_peak(intensity_polar, fargo_par)

    # Return peak relative to the center
    peak = theta[arg_phi]
    peak_offset = peak * (180.0 / np.pi) - 180.0

    # Store in mp_array
    peak_offsets[i] = peak_offset

def measure_peak_offset(args, threshold = 0.5):
    # Unpack
    i, frame = args

    intensity_polar = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number)
    peak_r_i, peak_phi_i = az.get_peak(intensity_polar, fargo_par)

    # Return peak relative to the center, where the edges are set by a threshold
    intensity_polar /= np.max(intensity_polar) # normalize

    intensity_polar_at_peak = intensity_polar[peak_r_i, :] # take azimuthal profile
    left_index = az.my_searchsorted(intensity_polar_at_peak, threshold)
    right_index = len(intensity_polar_at_peak) - az.my_searchsorted(intensity_polar_at_peak[::-1], threshold) - 1

    # Convert to theta
    left_theta = theta[left_index] * (180.0 / np.pi)
    right_theta = theta[right_index] * (180.0 / np.pi)
    center_theta = left_theta + (right_theta - left_theta) / 2.0
    
    peak_theta = theta[peak_phi_i] * (180.0 / np.pi)
    peak_offset = peak_theta - center_theta

    # Store in mp_array
    peak_offsets[i] = peak_offset


###############################################################################

##### PLOTTING #####

peak_offsets = mp_array("d", len(frame_range))

def make_plot(show = False):
    fig = plot.figure(figsize = (7, 5), dpi = dpi)
    ax = fig.add_subplot(111)

    # Get Data
    if num_cores == 1:
        if measure:
            peaks = np.array([measure_peak_offset((i, frame), threshold = threshold) for i, frame in enumerate(frame_range)])
        else:
            peaks = np.array([get_peak_offset((i, frame)) for i, frame in enumerate(frame_range)])
    else:
        # Pool
        if measure:
            offset_function = lambda f : measure_peak_offset(f, threshold = threshold)
        else:
            offset_function = get_peak_offset

        pool_args = [(i, frame) for (i, frame) in enumerate(frame_range)]

        p = Pool(num_cores)
        p.map(offset_function, pool_args)
        p.terminate()

    # Plot
    bins = np.linspace(-60, 60, 13) # Make this parameters
    data = np.array(peak_offsets)

    hist = plot.hist(data, bins = bins, c = 'r', histtype = 'step')
    hist_cum = plot.hist(data, bins = bins, c = 'b', histtype = 'step', cumulative = True)

    # Print
    if print_histogram:
        for i, (value, value_cum) in enumerate(zip(hist, hist_cum)):
            print "%f - %f: %d, %d, (%d, %d)" % (bins[i], bins[i + 1], value, value_cum, value / len(frame_range), value_cum / len(frame_range))

    # Save, Show, and Close
    frame_str = ""
    for i, frame_i in enumerate(args.frames):
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    if version is None:
        save_fn = "%s/intensityPeakHistogram_%s.png" % (save_directory, frame_str)
    else:
        save_fn = "%s/v%04d_intensityPeakHistogram_%s.png" % (save_directory, version, frame_str)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

    if save_data:
        save_name = "id%04d_b%02d_intensityPeaks.p" % (id_number, beam_size * planet_radius)
        pickle.dump(peaks, open(save_name, "wb"))

##### Make Plots! #####

make_plot(show = show)


