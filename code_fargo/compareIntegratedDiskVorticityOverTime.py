"""
compares excess masses with different taper times

Usage:
python compareExcessMassesOverTime.py
"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool
from multiprocessing import Array as mp_array

import math
import numpy as np
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot
from matplotlib import gridspec

from itertools import groupby

from pylab import rcParams # replace with rc ???
from pylab import fromfile

import util
from readTitle import readTitle

mass = 1
viscosity = -7
tapers = [10, 500, 1000, 2000]

##### PLOTTING #####

# Plot Parameters
linewidth = 4
fontsize = 14

my_dpi = 100
alpha = 0.5

max_frame = 0

# Plot Curves
for taper in tapers:
    frame_range = pickle.load(open("integrated_disk_vorticity_frames_taper%d.p" % taper, "rb"))
    vorticity_over_time = pickle.load(open("integrated_disk_vorticity_values_taper%d.p" % taper, "rb"))

    # make positive for log plot
    vorticity_over_time -= vorticity_over_time[0] # zero out start
    vorticity_over_time = np.abs(vorticity_over_time)

    # Curve
    plot.plot(frame_range, vorticity_over_time, linewidth = linewidth, label = r"$T_{growth}=$" + "%d" % taper)

    # Record Max Frame
    if frame_range[-1] > max_frame:
        max_frame = frame_range[-1]

# Annotate
title = r"$m_p = " + str(mass) + r" $ $M_J$, $\nu_{disk} = 10^{" + str(viscosity) + r"}$"
plot.xlabel("Number of Planet Orbits", fontsize = fontsize)
plot.ylabel("Integrated Disk Vorticity (+ Constant)", fontsize = fontsize)
plot.title(title, fontsize = fontsize + 2)

plot.legend(loc = "lower right")

# Axes
plot.xlim(10, max_frame)
#plot.ylim(0.0, 1.0)

plot.yscale("log")

# Save + Close
plot.savefig("integratedDiskVorticityOverTime_m%d_v%d.png" % (mass, abs(viscosity)))
plot.show()

plot.close()

