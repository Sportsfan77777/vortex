import numpy as np
from matplotlib import pyplot as plot

data = np.loadtxt("tqwk0.dat")

x = data[:,-1]
y1 = data[:,3] # Torque from Inner Disk 
y2 = data[:,4] # Torque from Outer Disk
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

