"""
convert dust density to cgs
--> can interpolate between sizes
--> can interpolate density to arbitary resolution

also creates associated files (radial, azimuthal, grain, temperature)
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

# Constants
G = 6.67 * 10**-8
mu = 2.34
mp = 1.67 * 10**-24
kb = 1.38 * 10**-16

# System Parameters
mass = 1.0 # (in solar masses)
radius = 5.0 # radius of planet (in AU)

mass_unit = mass * (1.988425 * 10**33) # (solar mass / g)
radius_unit = radius * (1.496 * 10**13) # (AU / cm)

density_unit = mass_unit / radius_unit**2 # unit conversion factor

# Grain Sizes
sizes = [0.01, 0.03, 0.1, 0.3, 1.0]
size_labels = ["hum", "hmm", "mm", "hcm", "cm"]

log_interpolated_sizes = np.linspace(np.log10(sizes[0]), np.log10(sizes[-1]), 100)
interpolated_sizes = np.power(10.0, log_interpolated_sizes) # to be used in size interpolation

# Output Resolution
new_num_rad = 300
new_num_theta = 400

# Save As
save_directory = "rt_input"
save_name = "double_interpolated"

######################################################################

# Input File
fn = sys.argv[1]

# Replace size with %s
for size in size_labels:
    fn = fn.replace(size, r"%s")

new_fn = "cgs_%s" % fn # output file

# Compile list of filenames
fns = []
for size in size_labels:
    fns.append(fn % size)

### Get FARGO Parameters ###
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = float(fargo_par["Nrad"])
num_theta = float(fargo_par["Nsec"])

rad = np.linspace(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]), num_rad + 1)[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta + 1)[:-1]

# Output Resolution
new_rad = np.linspace(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]), new_num_rad + 1)[:-1]
new_theta = np.linspace(0, 2 * np.pi, new_num_theta + 1)[:-1]

### Helper Functions ###
def find_argmax(density):
    # Returns Azimuthal Argmax
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    argmax = np.argmax(density_segment)
    arg_r, arg_phi = np.unravel_index(argmax, np.shape(density_segment))

    return arg_phi

def temperature(r):
    """ r (in cm) """
    # (mu * mp) (h**2 \omega **2) / (k_b)
    omega_sq = (G * mass_unit) / (r**3)
    scale_height_cgs = scale_height * (r)
    return (mu * mp) * (scale_height_cgs**2 * omega_sq) / (kb) 

# Get Data and Shift Vortex to Middle
middle = np.searchsorted(theta, np.pi)
density_arrays = {}

for (size_i, fn_i) in zip(size_labels, fns):
    tmp_density = fromfile(fn_i).reshape(num_rad, num_theta)
    location_i = find_argmax(tmp_density)

    shift_i = int(middle - location_i)
    density_arrays[size_i] = np.roll(tmp_density, shift_i, axis = 1)

# Convert Data and Interpolate to Arbitary Resolution
combination_array = np.zeros((new_num_rad * new_num_theta, len(sizes)))
for i, size_i in enumerate(size_labels):
    # Interpolate
    density_interpolation_function = sp_int.interp2d(rad, theta, density_arrays[size_i])
    interpolated_density = density_interpolation_function(new_rad, new_theta)

    # Convert to cgs
    combination_array[:, i] = (interpolated_density * density_unit).flatten('F') # F = column-major

# Interpolate to More Grain Sizes
size_interpolation_function = sp_int.interp1d(sizes, combination_array, axis = -1)
interpolated_combination_array = size_interpolation_function(interpolated_sizes)

# Combine (Interleave)
interleaved_array = interpolated_combination_array.flatten() # interleave to 1-d

# Bonus: Gather Temperature
temperatures = np.array([temperature(r * radius_unit) for r in new_rad])

# Save New Data
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

np.savetxt("%s/radial.dat" % save_directory, new_rad * radius_unit)
np.savetxt("%s/azimuthal.dat" % save_directory, new_theta)
np.savetxt("%s/grain.dat" % save_directory, interpolated_sizes)
np.savetxt("%s/temperature.dat" % save_directory, temperatures)
np.savetxt(save_directory + "/" + (new_fn % save_name), interleaved_array)

