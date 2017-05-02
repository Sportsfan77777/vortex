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
new_fn = "cgs_%s" % fn

# Replace size with %s
for size in sizes:
    fn = fn.replace(size, r"%s")

# Compile list of filenames
fns = []
new_fns = []
for size in sizes:
    fns.append(fn % size)

### Get FARGO Parameters ###
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = float(fargo_par["Nrad"])
num_theta = float(fargo_par["Nsec"])

# Data
density_arrays = {}
for (size_i, fn_i) in zip(sizes, fns):
    density_arrays[size_i] = fromfile(fn_i).reshape(num_rad, num_theta)

# Convert Data and Combine (Interleave) 
combination_array = np.zeros((num_rad * num_theta, len(sizes)))
for i, size_i in enumerate(sizes):
    combination_array[:, i] = (density_arrays[size_i] * density_unit).flatten('F') # F = column-major

combination_array = combination_array.flatten()

# Save New Data
np.save(new_fn, combination_array)






