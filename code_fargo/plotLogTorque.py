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

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
ks = 50.0 # Kernel Size
ks_small = ks / 3.0 # Smaller kernel to check the normal kernel

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt(torque_fn)
select = range(0, len(data[:,-1]), rate)
xs = (data[:,-1])[select] / (2 * np.pi) # Convert to num_orbits

# Load Planet Data for Analytic Comparison
planet_mass = 0.005 # In the future, load this from dictionary
viscosity = float(fargo_par["Viscosity"])
times = np.loadtxt(orbit_fn)[:, 0] # Note this is a different 'x' for plotting
radii = np.loadtxt(orbit_fn)[:, 2] # really 'a', not 'r'
def analytic_torque(rate = 50):
    """ returns time and torque at that time in [x, y] format """
    select = range(0, len(radii), rate)
    x = times[select] / (2 * np.pi)
    y = planet_mass * viscosity * np.power(np.array(radii[select]), -1.5)
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
    y4 = (np.sign(y2 - y1)) + 2 # Sign (inward is 3, outward is 1)
    y5 = analytic[1] # Analytic Reference

    s1 = np.abs(smooth(y1_base, ks)[select]) # Torque from Inner Disk (smoothed)
    s2 = np.abs(smooth(y2_base, ks)[select]) # Torque from Outer Disk (smoothed)
    s3 = np.abs(s2 - s1)
    s4 = (np.sign(s2 - s1)) + 2

    # Figure
    fig, (ax1, ax2) = plot.subplots(2, sharex = True)

    # Curves
    # Analytic
    plot.plot(analytic[0], analytic[1], c = "black", linewidth = linewidth)

    # Simulation
    ax1.plot(xs, y1, c = "r", alpha = alpha, linewidth = linewidth - 1)
    ax1.plot(xs, y2, c = "g", alpha = alpha, linewidth = linewidth - 1)
    ax1.plot(xs, y3, c = "b", alpha = alpha, linewidth = linewidth - 1)
    ax2.plot(xs, y4, c = "orange", alpha = alpha, linewidth = linewidth - 1)

    # Smoothed from Simulation
    ax1.plot(xs, s1, c = "r", label = "Inner", linewidth = linewidth)
    ax1.plot(xs, s2, c = "g", label = "Outer", linewidth = linewidth)
    ax1.plot(xs, s3, c = "b", label = "Inner + Outer", linewidth = linewidth)
    ax2.plot(xs, s4, c = "orange", linewidth = linewidth)

    ax1.plot([xs[0], xs[-1]], [0, 0], c = "black", linewidth = linewidth) # Zero Reference Line

    ax1.legend()

    # Annotate

    ax1.set_title(title, fontsize = fontsize + 2)
    ax2.set_xlabel("Number of Planet Orbits", fontsize = fontsize)
    ax1.set_ylabel("Torque", fontsize = fontsize)

    # Layout

    ax1.set_yscale('log')
    ax1.set_ylim(10**(-10), 10**(-3))

    ax2.set_ylim(0, 4) # Inward (3) or Outward (1)

    # Save and Close

    plot.savefig(save_fn, bbox_inches = 'tight')
    plot.show()

    plot.cla(fig)


### PLOTTING ###

make_plot(rla = False) # No Roche Lobe Avoidance
make_plot(rla = True) # Roche Lobe Tapering