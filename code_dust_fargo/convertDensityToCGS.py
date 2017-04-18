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

# Input File
fn = sys.argv[1]
new_fn = "cgs_%s" % fn

### Get FARGO Parameters ###
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = float(fargo_par["Nrad"])
num_theta = float(fargo_par["Nsec"])

# Data
density = fromfile(fn).reshape(num_rad, num_theta)

# Convert Data
density_cgs = density * density_unit

# Save New Data
np.save(new_fn, density_cgs)






