"""
Usage:
python generateSyntheticInput.py

Does
(1) Convert density distributions from code units to cgs
(2) Options to carve out inner cavity (default: yes) and scale dust sizes (default: no)
(3) Centers vortices (using cm, hcm, and mm peaks as a guide -- For T = 1000, center assumed to be fixed)
(4) Re-samples at a lower resolution
(5) Interpolates dust to continuous distribution of grain sizes (e.g. 100)
(6) Generates secondary input: radial.dat, azimuthal.dat, temperature.dat, grain.dat
(7) Outputs "gasddens%d.dat", a text file input for the synthetic image generator (access with "ln -s gasddens%d.dat" density.dat)
(8) Output "gassdens%d.p", a pickle file that can be used to look at the overall dust density distribution (such as to look at the gas-to-dust ratio over time)

Assumes directory structure with adjacent access to directories containing different grain sizes
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np
from scipy import interpolate as sp_int

import argparse
from pylab import fromfile

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

    # System Parameters
    parser.add_argument('-m', dest = "mass", type = float, default = 1.0,
                         help = 'mass of star in solar masses (default: 1.0)')
    parser.add_argument('-r', dest = "radius", type = float, default = 1.0,
                         help = 'radius of planet in AU (default: 20.0)')

    # Save Parameters
    parser.add_argument('-g', dest = "num_grains", type = int, default = 100,
                         help = 'number of interpolated grains (default: 100)')
    parser.add_argument('-s', dest = "new_res", nargs = 2, type = int, default = [300, 400],
                         help = 're-sample resolution (default: [300, 400])')
    parser.add_argument('--save', dest = "save_directory", default = ".",
                         help = 'save directory (default: ".")')

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

### Get Input Parameters ###

# Frames
if len(args.frames) == 1:
    frame_range = args.frames
elif len(args.frames) == 3:
    start = args.frames[0]; end = args.frames[1]; rate = args.frames[2]
    frame_range = range(start, end + 1, rate)

# Number of Cores 
num_cores = args.num_cores

# System Parameters
mass = args.mass # (in solar masses)
radius = args.radius # radius of planet (in AU)

# Save Parameters
new_num_rad = args.new_res[0]; new_num_theta = args.new_res[1]
save_directory = args.save_directory

### Miscellaneous ###
# Constants (G, \mu, m_p, k_b)
G = 6.67 * 10**-8; mu = 2.34; mp = 1.67 * 10**-24; kb = 1.38 * 10**-16

# Units
mass_unit = mass * (1.988425 * 10**33) # (solar mass / g)
radius_unit = radius * (1.496 * 10**13) # (AU / cm)

### Task Functions ###

def retrieve_density(frame, sizes):
    """ Step 0: Retrieve density """
    density = np.zeros((num_rad, num_theta, len(sizes)))

    for i, size_i in enumerate(sizes):
        fn_i = "../%s-size/gasddens%d.dat" % (size_i, frame)
        density[:, :, i] = fromfile(fn_i).reshape(num_rad, num_theta)

    return density

def convert_units(density):
    """ Step 1: convert density from code units to cgs units """
    density_unit = mass_unit / radius_unit**2 # unit conversion factor

    density *= density_unit
    return density

def polish(density, sizes, cavity_cutoff = 0.92, scale = 1):
    """ Step 2: get rid of inner cavity and scale dust densities to different grain size """
    # Cavity
    if cavity_cutoff is not None:
        cavity_cutoff_i = np.searchsorted(rad, cavity_cutoff)
        density[:cavity_cutoff_i, ] /= 100.0

    # Scale
    density *= scale
    sizes *= scale

    return density

def center_vortex(vortex_start = None, vortex_end = 2000):
    """ Step 3: center the vortex so that the peak is at 180 degrees """
    pass

def resample(new_num_rad = 300, new_num_theta = 400):
    """ Step 4: lower resolution (makes txt output smaller) """

    new_density = np.zeros((new_num_rad, new_num_theta, len(sizes)))

    new_rad = np.linspace(fargo_par["Rmin"], fargo_par["Rmax"], new_num_rad)
    new_theta = np.linspace(0, 2 * np.pi, new_num_theta)

    for i in len(sizes):
        interpolator = sp_int.interp2d(rad, theta, np.transpose(density[:, :, i])) # Careful: z is flattened!
        new_density[:, :, i] = (interpolator(new_rad, new_theta)).T # Note: result needs to be transposed
    
    return new_rad, new_theta, new_density

def interpolate_density(density, num_grains):
    """ Step 5: interpolate to more grain sizes """

    ### New Grain Sizes and Ranges ###
    log_interpolated_sizes = np.linspace(np.log10(min(sizes)), np.log10(max(sizes[0])), num_grains)
    interpolated_sizes = np.power(10.0, log_interpolated_sizes) # to be used in size interpolation

    # Grain Ranges (for power-law distribution)
    da = log_interpolated_sizes[1] - log_interpolated_sizes[0]
    log_interpolated_ranges = np.linspace(np.log10(min(sizes)) - 0.5 * da, np.log10(max(sizes)) + 0.5 * da, num_grains + 1)
    interpolated_ranges = np.power(10.0, log_interpolated_ranges)

    ### Interpolate ### 
    size_interpolator = sp_int.interp1d(sizes, density, axis = -1)
    unweighted_interpolated_density = size_interpolation_function(interpolated_sizes)

    # Scale to a Power Law Distribution
    def integrate_size_distribution(x, y):
        """ Determine weight of each dust size bin """
        def f(z):
            return z ** (surface_density_power + 1.0) # +1 for integration
        return f(y) - f(x)

    range_a = interpolated_ranges[:-1]
    range_b = interpolated_ranges[1:]
    size_weights = [integrate_size_distribution(x, y) for (x, y) in zip(range_a, range_b)]
    power_law = size_weights / np.sum(size_weights)

    interpolated_density = power_law[None, :] * unweighted_interpolated_density

    return interpolated_density

def generate_secondary_files(rad, theta, sizes):
    """ Step 6: write other *.dat files """

    def temperature(R):
        """ r """
        R *= radius_unit # convert to cgs

        # (mu * mp) (h**2 \omega **2) / (k_b)
        omega_sq = (R * mass_unit) / (R**3)
        scale_height_cgs = scale_height * (R)
        return (mu * mp) * (scale_height_cgs**2 * omega_sq) / (kb) 
        
    temperatures = np.array([temperature(r) for r in rad])

    # Save Files
    np.savetxt("%s/radial.dat" % save_directory, rad * radius_unit)
    np.savetxt("%s/azimuthal.dat" % save_directory, theta)
    np.savetxt("%s/grain.dat" % save_directory, sizes)
    np.savetxt("%s/temperature.dat" % save_directory, temperatures)

def output_density_txt(density, frame):
    """ Step 7: output txt file (combine with step 8?) """
    interleaved_density = density.flatten('F') # interleave to 1-d

    fn = "gasddens%d.dat" % frame
    np.savetxt(save_directory + "/" + (new_fn % save_name), interleaved_density)

def output_density_pickle(density, frame):
    """ Step 8: output pickle file """
    fn = "gasddens%d.p" % frame
    pickle.dump(density, open(fn, 'w'))

def full_procedure(frame):
    """ Every Step """
    density = retrieve_density(frame, sizes)

    density = convert_units(density)
    density = polish(density)
    density = center_vortex(density)
    new_rad, new_theta, density = resample(new_num_rad = new_num_rad, new_num_theta = new_num_theta)
    density = interpolate_density(density, num_grains)
    output_density_txt(density, frame)
    output_density_pickle(density, frame)

    if frame == frame_range[0]:
        generate_secondary_files(new_rad, new_theta, sizes)


### Generate Synthetic Input ###

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


