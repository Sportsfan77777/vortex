"""
plots torque over time

 *** uses a smoothing function ***
"""

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc

from scipy import signal as sig
from scipy.ndimage import filters as ff

# Smoothing Function
smooth = lambda array, kernel_size : ff.gaussian_filter(array, kernel_size) # smoothing filter
ks = 50.0 # Kernel Size
ks_small = ks / 4.0 # Smaller kernel to check the normal kernel

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt("tqwk0.dat")
select = range(0, len(data[:,-1]), rate)
x = (data[:,-1])[select]

# Plot Parameters
alpha = 0.2 # for non-smoothed curves
fontsize = 14
linewidth = 2

def make_plot(rla = True):
    if not rla:
        # No Roche Lobe Avoidance
        y1_base = data[:,1]
        y2_base = data[:,2]
        title = "Torque! (with no roche lobe avoidance)"
        save_fn = "TorqueEvolution_with_no_rla.png"
    else:
        # Roche Lobe Tapering
        y1_base = data[:,3]
        y2_base = data[:,4]
        title = "Torque! (with roche lobe tapering)"
        save_fn = "TorqueEvolution_with_rla.png"

    # Data

    y1 = smooth(y1_base, ks_small)[select] # Torque from Inner Disk 
    y2 = smooth(y2_base, ks_small)[select] # Torque from Outer Disk
    y3 = y1 + y2

    s1 = smooth(y1_base, ks)[select] # Torque from Inner Disk (smoothed)
    s2 = smooth(y2_base, ks)[select] # Torque from Outer Disk (smoothed)
    s3 = s1 + s2

    # Curves

    plot.plot(x, y1, c = "r", alpha = alpha, linewidth = linewidth - 1)
    plot.plot(x, y2, c = "g", alpha = alpha, linewidth = linewidth - 1)
    plot.plot(x, y3, c = "b", alpha = alpha, linewidth = linewidth - 1)

    plot.plot(x, s1, c = "r", label = "Inner", linewidth = linewidth)
    plot.plot(x, s2, c = "g", label = "Outer", linewidth = linewidth)
    plot.plot(x, s3, c = "b", label = "Inner + Outer", linewidth = linewidth)

    plot.plot([x[0], x[-1]], [0, 0], c = "black", linewidth = linewidth) # Zero Reference Line

    plot.legend()

    # Annotate

    plot.title(title, fontsize = fontsize + 2)
    plot.xlabel("Timestep", fontsize = fontsize)
    plot.ylabel("Torque", fontsize = fontsize)

    # Save and Close

    plot.savefig(save_fn, bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot(rla = True) # Roche Lobe Tapering
make_plot(rla = False) # No Roche Lobe Avoidance