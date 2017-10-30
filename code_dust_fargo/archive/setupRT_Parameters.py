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

### Get FARGO Parameters ###
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

scale_height = float(fargo_par["AspectRatio"])

num_rad = float(fargo_par["Nrad"])
num_theta = float(fargo_par["Nsec"])

r_in = float(fargo_par["Rin"])
r_out = float(fargo_par["Rout"])

# Helper Functions
def temperature(r):
    """ r (in cm) """
    # (mu * mp) (h**2 \omega **2) / (k_b)
    omega_sq = (G * mass_unit) / (r**3)
    scale_height_cgs = scale_height * (r)
    return (mu * mp) * (scale_height_cgs**2 * omega_sq) / (kb) 

# Parameters
radii = np.linspace(r_in, r_out, num_rad + 1) * (radius_unit)
angles = np.linspace(0, 2 * np.pi, num_theta)
temperatures = np.array([temperature(r) for r in radii]) 

# Save New Data
np.savetxt("radial.dat", radii[:-1])
np.savetxt("azimuthal.dat", angles)
np.savetxt("temperature.dat", temperatures[:-1])

