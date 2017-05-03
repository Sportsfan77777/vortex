"""
convert dust density to cgs
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np

from pylab import fromfile

# System Parameters
mass = 1.0 # (in solar masses)
radius = 5.0 # radius of planet (in AU)

mass_unit = mass * (1.988425 * 10**33) # (solar mass / g)
radius_unit = radius * (1.496 * 10**13) # (AU / cm)

density_unit = mass_unit / radius_unit**3 # unit conversion factor

# Grain Sizes
sizes = ["cm", "hcm", "mm", "hmm", "hum"]

######################################################################

# Input File
fn = sys.argv[1]

# Replace size with %s
for size in sizes:
    fn = fn.replace(size, r"%s")

new_fn = "cgs_%s" % fn # output file

# Compile list of filenames
fns = []
for size in sizes:
    fns.append(fn % size)

### Get FARGO Parameters ###
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = float(fargo_par["Nrad"])
num_theta = float(fargo_par["Nsec"])

rad = np.linspace(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]), num_rad + 1)
theta = np.linspace(0, 2 * np.pi, num_theta + 1)

### Helper Function ###
def find_argmax(density):
    # Returns Azimuthal Argmax
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    argmax = np.argmax(density_segment)
    arg_r, arg_phi = np.unravel_index(argmax, np.shape(density_segment))

    return arg_phi

# Get Data and Shift Vortex to Middle
middle = np.searchsorted(theta, np.pi)
density_arrays = {}

for (size_i, fn_i) in zip(sizes, fns):
    tmp_density = fromfile(fn_i).reshape(num_rad, num_theta)
    location_i = find_argmax(tmp_density)

    shift_i = int((middle - location_i) * (num_theta))
    density_arrays[size_i] = np.roll(tmp_density, shift_i, axis = 1)

# Convert Data and Combine (Interleave) 
combination_array = np.zeros((num_rad * num_theta, len(sizes)))
for i, size_i in enumerate(sizes):
    combination_array[:, i] = (density_arrays[size_i] * density_unit).flatten('F') # F = column-major

combination_array = combination_array.flatten() # interleave to 1-d

# Save New Data
np.savetxt(new_fn % "all", combination_array)

