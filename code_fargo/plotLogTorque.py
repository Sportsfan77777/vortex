"""
plots torque over time

 *** uses a smoothing function ***
"""

import sys
import os
import subprocess

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

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
ks = 50.0 # Kernel Size
ks_small = ks / 3.0 # Smaller kernel to check the normal kernel

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt("tqwk0.dat")
select = range(0, len(data[:,-1]), rate)
x = (data[:,-1])[select] / (2 * np.pi) # Convert to num_orbits

# Load Planet Data for Analytic Comparison
planet_mass = 0.005 # In the future, load this from dictionary
viscosity = float(fargo_par["Viscosity"])
times = np.loadtxt("orbit.dat")[0] # Note this is a different 'x' for plotting
radius = np.loadtxt("orbit.dat")[2] # really 'a', not 'r'
def analytic_torque(rate = 10):
    """ returns time and torque at that time in [x, y] format """
    select = range(0, len(radius), rate)
    x = times[select] / (2 * np.pi)
    y = planet_mass * viscosity * np.power(np.array(radius[select]), -1.5)
    # y = torque = J / t = (M R^2 \Omega) / (R^2 / \nu)
    # y = M * \nu * \Omega
    return [x, y]

# Plot Parameters
alpha = 0.2 # for non-smoothed curves
fontsize = 14
linewidth = 2

def make_plot(rla = True):
    if not rla:
        # No Roche Lobe Avoidance
        y1_base = data[:,1]
        y2_base = data[:,2]
        title = "Torque (0th): %s x %s" % (fargo_par["Nrad"], fargo_par["Nsec"])
        save_fn = "TorqueZerothOrder.png"
    else:
        # Roche Lobe Tapering
        y1_base = data[:,3]
        y2_base = data[:,4]
        title = "Torque (1st): %s x %s" % (fargo_par["Nrad"], fargo_par["Nsec"])
        save_fn = "TorqueFirstOrder.png"

    # Data
    analytic = analytic_torque()

    y1 = np.abs(smooth(y1_base, ks_small)[select]) # Torque from Inner Disk 
    y2 = np.abs(smooth(y2_base, ks_small)[select]) # Torque from Outer Disk
    y3 = np.abs(y2 - y1) # Difference
    y4 = (np.sign(y2 - y1) * 10**(-9)) + 3 * 10**(-9) # Sign
    y5 = analytic[1] # Analytic Reference

    s1 = np.abs(smooth(y1_base, ks)[select]) # Torque from Inner Disk (smoothed)
    s2 = np.abs(smooth(y2_base, ks)[select]) # Torque from Outer Disk (smoothed)
    s3 = np.abs(s2 - s1)
    s4 = (np.sign(s2 - s1) * 10**(-9)) + 3 * 10**(-9)

    # Curves
    # Analytic
    plot.plot(analytic[0], analytic[1], c = "black", linewidth = linewidth - 1)

    # Simulation
    plot.plot(x, y1, c = "r", alpha = alpha, linewidth = linewidth - 1)
    plot.plot(x, y2, c = "g", alpha = alpha, linewidth = linewidth - 1)
    plot.plot(x, y3, c = "b", alpha = alpha, linewidth = linewidth - 1)
    plot.plot(x, y4, c = "orange", alpha = alpha, linewidth = linewidth - 1)

    # Smoothed from Simulation
    plot.plot(x, s1, c = "r", label = "Inner", linewidth = linewidth)
    plot.plot(x, s2, c = "g", label = "Outer", linewidth = linewidth)
    plot.plot(x, s3, c = "b", label = "Inner + Outer", linewidth = linewidth)
    plot.plot(x, s4, c = "orange", linewidth = linewidth)


    plot.plot([x[0], x[-1]], [0, 0], c = "black", linewidth = linewidth) # Zero Reference Line

    plot.legend()

    # Annotate

    plot.title(title, fontsize = fontsize + 2)
    plot.xlabel("Timestep", fontsize = fontsize)
    plot.ylabel("Torque", fontsize = fontsize)

    # Layout

    plot.yscale('log')
    plot.ylim(10**(-10), 10**(-4))

    # Save and Close

    plot.savefig(save_fn, bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot(rla = True) # Roche Lobe Tapering
make_plot(rla = False) # No Roche Lobe Avoidance