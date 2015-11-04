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
smooth = lambda array : ff.gaussian_filter(array, 750) # choose smoothing filter

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt("tqwk0.dat")
select = range(0, len(data[:,-1]), rate)
x = (data[:,-1])[select]

# Plot Parameters
alpha = 0.2 # for non-smoothed curves
fontsize = 14

############# NO ROCHE LOBE AVOIDANCE #############

# Data

y1 = smooth(data[:,1])[select] # Torque from Inner Disk 
y2 = smooth(data[:,2])[select] # Torque from Outer Disk
y3 = y1 + y2

s1 = smooth(data[:,1])[select] # Torque from Inner Disk 
s2 = smooth(data[:,2])[select] # Torque from Outer Disk
s3 = s1 + s2

# Curves

plot.plot(x, y1, c = "r", alpha = alpha)
plot.plot(x, y2, c = "g", alpha = alpha)
plot.plot(x, y3, c = "b", alpha = alpha)

plot.plot(x, s1, c = "r", label = "Inner")
plot.plot(x, s2, c = "g", label = "Outer")
plot.plot(x, s3, c = "b", label = "Inner + Outer")

plot.legend()

# Annotate

plot.title("Torque! (with no roche lobe avoidance)", fontsize = fontsize + 2)
plot.xlabel("Time", fontsize = fontsize)
plot.ylabel("Torque", fontsize = fontsize)

# Save and close

plot.savefig("TorqueEvolution_no_rla.png", bbox_inches = 'tight')
plot.show()

plot.cla()

############# ROCHE LOBE TAPERING #############

# Data

y1 = (data[:,3])[select] # Torque from Inner Disk 
y2 = (data[:,4])[select] # Torque from Outer Disk
y3 = y1 + y2

s1 = smooth(data[:,3])[select] # Torque from Inner Disk 
s2 = smooth(data[:,4])[select] # Torque from Outer Disk
s3 = s1 + s2

# Curves

plot.plot(x, y1, c = "r", alpha = alpha)
plot.plot(x, y2, c = "g", alpha = alpha)
plot.plot(x, y3, c = "b", alpha = alpha)

plot.plot(x, s1, c = "r", label = "Inner Disk")
plot.plot(x, s2, c = "g", label = "Outer Disk")
plot.plot(x, s3, c = "b", label = "Inner + Outer")

plot.legend()

# Annotate

plot.title("Torque! (with roche lobe tapering)", fontsize = fontsize + 2)
plot.xlabel("Time", fontsize = fontsize)
plot.ylabel("Torque (units???)", fontsize = fontsize)

# Save

plot.savefig("TorqueEvolution_with_rla.png", bbox_inches = 'tight')
plot.show()
