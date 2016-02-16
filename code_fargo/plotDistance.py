"""
plots torque over time

 *** uses a smoothing function ***
"""

import sys
import os
import subprocess
import pickle

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc

from scipy import signal as sig
from scipy.ndimage import filters as ff

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt("orbit0.dat")
select = range(0, len(data[:,-1]), rate)
xs = (data[:,0])[select] / (2 * np.pi) # Convert to num_orbits
sm_axes = (data[:,2])[select] # Planet Semi-Major Axis
eccs = (data[:,1])[select] # Planet Eccentricity

# Calculate Analytic Rate
viscosity = float(fargo_par["Viscosity"])
delta_t = data[1,0] - data[0,0]
def next_a(previous_a, delta_t):
    timescale = previous_a**2 / viscosity
    return previous_a - delta_t / timescale

ys_analytic = [sm_axes[0]]
for i, x in enumerate(xs[:-1]):
    a = next_a(ys_analytic[-1], delta_t) # Iterate through time
    ys_analytic.append(a)

# Plot Parameters
alpha = 0.2 # for non-smoothed curves
fontsize = 14
linewidth = 3

def make_plot():
    # Curves
    ys = sm_axes
    qs = np.array([a * (1 - e) for (a, e) in zip(ys, eccs)])

    # Data
    plot.plot(xs, ys, c = "blue", linewidth = linewidth, label = "a")
    plot.plot(xs, qs, c = "purple", alpha = alpha, linewidth = linewidth, label = "q")

    # Analytic
    plot.plot(xs, ys_analytic, c = "black", linewidth = linewidth, label = "ideal")

    # Annotate
    plot.title("Distance Over Time", fontsize = fontsize + 2)
    plot.xlabel("Timestep", fontsize = fontsize)
    plot.ylabel("Planet Distance ('a' or 'q')", fontsize = fontsize)

    plot.legend()

    # Limits
    mins = [min(ys), min(qs), min(ys_analytic)]
    maxes = [max(ys), max(qs), max(ys_analytic)]
    plot.ylim(min(mins)- 0.002, max(maxes) + 0.025)

    # Save and Close
    plot.savefig("planetDistance.png", bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot()