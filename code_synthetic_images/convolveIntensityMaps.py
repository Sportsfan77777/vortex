"""
convolves intensity map to make it look like an alma image

python convolveIntensityMap.py frame
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

from scipy import signal
from scipy.interpolate import interp1d as interpolate
from scipy.ndimage import map_coordinates

import matplotlib
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
import square as sq

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot convolved intensity maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = None,
                         help = 'save directory (default: beam%03d % beam_size in AU)')

    # Convolution Parameters
    parser.add_argument('-b', dest = "beam_size", type = float, default = 0.5,
                         help = 'beam size (in planet radii) (default: 0.5)')

    # Identification
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get ID%04d Parameters ###
fn = "id%04d_par.p" % args.id_number
fargo_par = pickle.load(open(fn, "rb"))

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

### Get Synthetic Image Parameters ###
synthetic_par = np.loadtxt("../parameters.dat")
wavelength = int(float(synthetic_par[2]))
distance = float(synthetic_par[4])

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Convolution Parameters
beam_size = args.beam_size
beam_radius = 0.5 * beam_size # the sigma of the gaussian, not the beam diameter
beam_diameter = beam_size * fargo_par["Radius"]

# Files
save_directory = args.save_directory
if save_directory is None:
    save_directory = "beam%03d" % (beam_diameter)
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Identification
id_number = args.id_number
version = args.version

###############################################################################

### Add new parameters to dictionary ###
fargo_par["Beam"] = beam_size
fargo_par["Wavelength"] = wavelength
fargo_par["Distance"] = distance

###############################################################################

### Helper Functions ###

def read_data(frame):
    """ Step 0: read intensity data """
    intensity = util.read_data(frame, 'intensity', fargo_par, id_number = id_number)
    return intensity

def clear_inner_disk(intensity):
    """ Step 1: get rid of inner disk (r < outer_limit) """
    filtered_intensity = np.copy(intensity)

    outer_limit = np.searchsorted(rad, 0.75) # make outer limit a parameter
    filtered_intensity[:outer_limit] = 0

    return filtered_intensity

######## Add zeroed-out inner disk!!!!!! #########

def polar_to_cartesian(intensity, order = 3):
    """ Step 2: convert to cartesian """
    xs, ys, _, _, cartesian_intensity = sq.polar_to_cartesian(intensity, rad, theta)
    return xs, ys, cartesian_intensity

def convolve_intensity(intensity):
    """ Step 3: Convolve cartesian intensity with 2-D Gaussian """
    # Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.fftconvolve.html#scipy.signal.fftconvolve

    # Determine Gaussian Parameters
    dr = rad[1] - rad[0]
    sigma = int(beam_radius / dr)
    window = signal.gaussian(5 * sigma, std = sigma)

    # Construct 2-D Gaussian Kernel
    normed_window = window / np.sum(window)
    kernel = np.outer(normed_window, normed_window)

    convolved_intensity = signal.fftconvolve(intensity, kernel, mode = 'same')
    return convolved_intensity

def divide_by_beam(intensity):
    """ Step 4: divide by beam """
    # beam size = (pi / 4 ln2) * (theta)^2
    beam_angle = beam_diameter / distance
    beam = np.pi * (beam_angle)**2 / (4 * np.log(2))

    return intensity / beam

def save_data(intensity_cart, frame, xs, ys, order = 3):
    """ Step 5: save in cartesian and polar coordinates """

    ################### Save cartesian output ###################
    # Save in pickle
    if version is None:
        save_fn = "%s/id%04d_intensityCartGrid%04d.p" % (save_directory, id_number, frame)
    else:
        save_fn = "%s/v%04d_id%04d_intensityCartGrid%04d.p" % (save_directory, version, id_number, frame)
    pickle.dump(intensity_cart, open(save_fn, 'wb'))

    ##################### Save polar output #####################
    # Set up rt-grid
    old_rs = rad
    old_thetas = theta

    rs_grid, thetas_grid = np.meshgrid(old_rs, old_thetas)

    # Interpolate xy-grid
    interpolated_xs = interpolate(xs, np.arange(len(xs)), bounds_error = False)
    interpolated_ys = interpolate(ys, np.arange(len(ys)), bounds_error = False)

    # Match up rt-grid with xy-grid
    new_xs = rs_grid * np.cos(thetas_grid)
    new_ys = rs_grid * np.sin(thetas_grid)

    new_interpolated_xs = interpolated_xs(new_xs.ravel())
    new_interpolated_ys = interpolated_ys(new_ys.ravel())

    # Convert (Interpolate) and Transpose
    polar_data = map_coordinates(intensity_cart, np.array([new_interpolated_xs, new_interpolated_ys]), order = order).reshape(new_xs.shape)
    polar_data = polar_data.T

    # Rotate 90 degrees! (artifact of method?)
    quarter_turn = -np.searchsorted(theta, np.pi / 2.0)
    polar_data = np.roll(polar_data, quarter_turn, axis = -1)

    # Save in pickle
    if version is None:
        save_fn = "%s/id%04d_intensityMap%04d.p" % (save_directory, id_number, frame)
    else:
        save_fn = "%s/v%04d_id%04d_intensityMap%04d.p" % (save_directory, version, id_number, frame)
    pickle.dump(polar_data, open(save_fn, 'wb'))

    # Save fargo_par
    dict_name = "%s/id%04d_par.p" % (save_directory, id_number)
    pickle.dump(fargo_par, open(dict_name, "wb"))

##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

def full_procedure(frame):
    """ Every Step """

    intensity = read_data(frame)
    intensity = clear_inner_disk(intensity) # not neccessary if inner disk was never sampled
    xs, ys, intensity_cartesian = polar_to_cartesian(intensity)

    intensity_cartesian = convolve_intensity(intensity_cartesian)
    intensity_cartesian = divide_by_beam(intensity_cartesian)
    save_data(intensity_cartesian, frame, xs, ys)


##### Make Plots! #####

# Iterate through frames

if len(frame_range) == 1:
    full_procedure(frame_range[0])
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(full_procedure, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            full_procedure(frame)



