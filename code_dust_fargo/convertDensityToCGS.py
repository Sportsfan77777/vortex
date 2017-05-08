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

density_unit = mass_unit / radius_unit**3 # unit conversion factor

# Output Resolution
new_num_rad = 300
new_num_theta = 400

# Save As
save_directory = "rt_input"

######################################################################

# Input File
fn = sys.argv[1]
new_fn = "cgs_%s" % fn

### Get FARGO Parameters ###
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

scale_height = float(fargo_par["AspectRatio"])

num_rad = float(fargo_par["Nrad"])
num_theta = float(fargo_par["Nsec"])

rad = np.linspace(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]), num_rad + 1)[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta + 1)[:-1]

# Output Resolution
new_rad = np.linspace(float(fargo_par["Rmin"]), float(fargo_par["Rmax"]), new_num_rad + 1)[:-1]
new_theta = np.linspace(0, 2 * np.pi, new_num_theta + 1)[:-1]

# Helper Functions
def temperature(r):
    """ r (in cm) """
    # (mu * mp) (h**2 \omega **2) / (k_b)
    omega_sq = (G * mass_unit) / (r**3)
    scale_height_cgs = scale_height * (r)
    return (mu * mp) * (scale_height_cgs**2 * omega_sq) / (kb) 

### Data ###
density = fromfile(fn).reshape(num_rad, num_theta)

# Interpolate to Arbitary Resolution
density_interpolation_function = sp_int.interp2d(rad, theta, density)
interpolated_density = density_interpolation_function(new_rad, new_theta)

# Convert Data
density_cgs = interpolated_density * density_unit

# Bonus: Gather Temperature
temperatures = np.array([temperature(r) for r in new_rad])

# Save New Data
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

np.savetxt("%s/radial.dat" % save_directory, new_rad)
np.savetxt("%s/azimuthal.dat" % save_directory, new_theta)
np.savetxt("%s/grain.dat" % save_directory, interpolated_sizes)
np.savetxt("%s/temperature.dat" % save_directory, temperatures)
np.save(save_directory + new_fn, density_cgs)






