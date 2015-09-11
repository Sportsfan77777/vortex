"""
plots torque over time

 *** uses a smoothing function ***
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot

from scipy import signal as sig
from scipy.ndimage import filters as ff
smooth = lambda array : ff.gaussian_filter(array, 750) # choose smoothing filter
ks = 1000 # choose kernel size

data = np.loadtxt("tqwk0.dat")

select = range(0, len(data[:,-1]), 4000)
x = (data[:,-1])[select]

x0 = data[:,-1]
x1 = np.roll(x0, 1)
d = x1 - x0
d[d > 0] = 0
d[d < 0] = 1

for i, q in enumerate(x0[125000:125500]):
	print i, q

print np.sum(d), len(d)

############# NO ROCHE LOBE AVOIDANCE #############

y1 = smooth(data[:,1])[select] # Torque from Inner Disk 
y2 = smooth(data[:,2])[select] # Torque from Outer Disk
y3 = y1 + y2

plot.plot(x, y1, c = "r", label = "Inner Disk")
plot.plot(x, y2, c = "g", label = "Outer Disk")
plot.plot(x, y3, c = "b", label = "Inner + Outer")

plot.legend()

fontsize = 14

plot.title("Torque! (with no roche lobe avoidance)", fontsize = fontsize + 2)
plot.xlabel("Time", fontsize = fontsize)
plot.ylabel("Torque", fontsize = fontsize)

plot.savefig("TorqueEvolution_no_roche_lobe_avoidance.png", bbox_inches = 'tight')
plot.show()

plot.cla()

############# ROCHE LOBE TAPERING #############

y1 = smooth(data[:,3])[select] # Torque from Inner Disk 
y2 = smooth(data[:,4])[select] # Torque from Outer Disk
y3 = y1 + y2

plot.plot(x, y1, c = "r", label = "Inner Disk")
plot.plot(x, y2, c = "g", label = "Outer Disk")
plot.plot(x, y3, c = "b", label = "Inner + Outer")

plot.legend()

fontsize = 14

plot.title("Torque! (with roche lobe tapering)", fontsize = fontsize + 2)
plot.xlabel("Time", fontsize = fontsize)
plot.ylabel("Torque (units???)", fontsize = fontsize)

plot.savefig("TorqueEvolution_with_roche_lobe_tapering.png", bbox_inches = 'tight')
plot.show()
