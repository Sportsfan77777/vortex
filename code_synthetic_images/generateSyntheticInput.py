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

from pylab import fromfile

directories = ["cm", "hcm", "mm", "hmm", "hum", "um"]
sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001])

# System Parameters
mass = 1.0 # (in solar masses)
radius = 20.0 # radius of planet (in AU)

# Constants (G, \mu, m_p, k_b)
G = 6.67 * 10**-8; mu = 2.34; mp = 1.67 * 10**-24; kb = 1.38 * 10**-16

# Units
mass_unit = mass * (1.988425 * 10**33) # (solar mass / g)
radius_unit = radius * (1.496 * 10**13) # (AU / cm)

def retrieve_density(sizes):
    """ Step 0: Retrieve density """
    density = np.zeros((num_rad * num_theta, len(sizes)))

    for i, (size_i, fn_i) in enumerate(zip(size_labels, fns)):
        density[:, :, i] = fromfile(fn_i).reshape(num_rad, num_theta)

    return density

def convert_units(density):
    """ Step 1: convert density from code units to cgs units """
    density_unit = mass_unit / radius_unit**2 # unit conversion factor

    density *= density_unit

def polish(density, sizes, cavity_cutoff = 0.92, scale = 1):
    """ Step 2: get rid of inner cavity and scale dust densities to different grain size """
    # Cavity
    if cavity_cutoff is not None:
        cavity_cutoff_i = np.searchsorted(rad, cavity_cutoff)
        density[:cavity_cutoff_i, ] /= 100.0

    # Scale
    density *= scale
    sizes *= scale

def center_vortex(vortex_start = None, vortex_end = 2000):
    """ Step 3: center the vortex so that the peak is at 180 degrees """
    pass

def resample(new_num_rad = 300, new_num_theta = 400):
    """ Step 4: lower resolution (makes txt output smaller) """

    new_rad = np.linspace(fargo_par["Rmin"], fargo_par["Rmax"], new_num_rad)
    new_theta = np.linspace(0, 2 * np.pi, new_num_theta)

    # Note: how to handle density arrays??? 6 separate, or 1 with an extra dimension
    # The input to interp2d must be 2-D!
    interpolator = sp_int.interp2d(rad, theta, np.transpose(density_arrays[size_i])) # Careful: z is flattened!
    new_density = (interpolator(new_rad, new_theta)).T # Note: result needs to be transposed
    
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

def output_density_txt(density):
    """ Step 7: output txt file (combine with step 8?) """
    interleaved_density = power_law_array.flatten('F') # interleave to 1-d
    np.savetxt(save_directory + "/" + (new_fn % save_name), interleaved_density)

def output_density_pickle(density):
    """ Step 8: output pickle file """

    #pickle.dump(density, open(, ))
    pass

### Generate Synthetic Input ###

density = retrieve_density()

convert_units(density)
polish(density)
center_vortex()
resample()
interpolate_density()
generate_secondary_files()
output_density_txt()
output_density_pickle()
