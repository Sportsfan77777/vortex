"""
scatter plot of azimuthal extents

Usage:
python scatterExtents.py
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

# Dictionaries
cases = [1, 5, 10, 50]
tapers = {}
extents = {}

# M = 1, v = 10^-6
case = 1
tapers[case] = [10, 250, 500]
extents[case] = [180, 240, 240]

# M = 1, v = 10^-7
case = 10
tapers[case] = [10, 500, 1000, 2000]
extents[case] = [120, 240, 240, 240]

# M = 5, v = 10^-6
case = 5
tapers[case] = [2, 100, 200]
extents[case] = [90, 120, 180]

# M = 5, v = 10^-7
case = 50
tapers[case] = [2, 200, 400, 800]
extents[case] = [120, 120, 120, 180]

#### PLOTTING ####
fontsize = 16

for case in cases:
	tapers_i = tapers[case]
	extents_i = extents[case]

	for i, (taper, extent) in enumerate(zip(tapers_i, extents_i)):
		if extent > 180:
			extent_type = "x"
		elif extent == 180:
			extent_type = "_"
		else:
			extent_type = "o"

		plot.scatter([case], [taper], marker = extent_type, s = 100)

# Annotate
plot.xlabel(r"$q / \nu$ $\times$ $[10^{6}]$", fontsize = fontsize)
plot.ylabel(r"$T_{growth}$ to Jupiter-size", fontsize = fontsize)
plot.title("Azimuthal Extents", fontsize = fontsize + 2)

# Axes
plot.xlim(0.5, 100)
#plot.ylim(1, 3000)

plot.xscale("log")
#plot.yscale("log")

# Save, Show, and Close
plot.savefig("scattered_extents.png")
plot.show()
plot.cla()


