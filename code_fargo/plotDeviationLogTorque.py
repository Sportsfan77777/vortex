"""
plots standard deviation of outer disk torque over time
"""

import sys
import os
import subprocess
import pickle

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plot
from matplotlib import gridspec
from matplotlib import rcParams as rc

from scipy import signal as sig
from scipy.ndimage import filters as ff

from readTitle import readTitle

## Set file names ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    torque_fn = "tqwk1.dat"
    orbit_fn = "orbit1.dat"
else:
    # fargo
    torque_fn = "tqwk0.dat"
    orbit_fn = "orbit0.dat"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

planet_mass = float(fargo_par["PlanetMass"])

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
ks = 50.0 # Kernel Size
ks_small = 1.0 # ks / 3.0 # Smaller kernel to check the normal kernel

# Load Data (and choose subset) = x-axis
rate = 5 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt(torque_fn)
select = range(0, len(data[:, -1]), rate)
xs = (data[:,-1])[select] / (2 * np.pi) # Convert to num_orbits

### Data ###
## For Comparison (Smoothed)
ks = 50
smooth_torque = smooth(data[:, 4], ks)
median_torque = np.median(smooth_torque)

## For Deviations (Sampled)
outer_disk_torque = planet_mass * (data[:, 4])[select] # Column 4 is the Outer Disk Torque with Roche Lobe Tapering
log_outer_disk_torque = np.log(outer_disk_torque) / np.log(10.0)

half_width = 5
start = half_width
end = len(outer_disk_torque) - half_width

torque_deviations = []
for center in range(start, end):
    deviation = np.std(outer_disk_torque[center - half_width : center + half_width])
    torque_deviations.append(deviation)

#torque_deviations = [10**(d) for d in torque_deviations]

# Plot Parameters
alpha = 0.2 # for non-smoothed curves
fontsize = 14
linewidth = 2

def make_plot(rla = True):
    # Figure
    fig = plot.figure(figsize = (8, 6)) 

    # Curves
    plot.plot(xs, outer_disk_torque, color = "blue", linewidth = linewidth, alpha = alpha)
    plot.plot(xs[start : end], torque_deviations, color = "red", linewidth = linewidth)
    plot.plot(xs[start : end], (torque_deviations / smooth_torque) * median_torque, color = "green", linewidth = linewidth)
    
    # Limits
    plot.xlim(0, xs[-1])

    plot.yscale('log')
    plot.ylim(10**(-11), 10**(-5))

    # Annotate
    this_title = readTitle()
    plot.title(this_title, fontsize = fontsize + 2)
    plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
    plot.ylabel("Torque", fontsize = fontsize)

    # Save and Close

    plot.savefig("torqueDeviation.png", bbox_inches = 'tight')
    plot.show()

    plot.close(fig)


### PLOTTING ###

make_plot()