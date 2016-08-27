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
cases = [1, 10, 25, 250, 2.01, 4.01, 10.01, 25.01, 75.01]
tapers = {}
extents = {}

# M = 1, v = 10^-6
case = cases[0]
tapers[case] = [10, 50, 100, 150, 250, 500]
extents[case] = [180, 180, 180, 240, 240, 240]

# M = 1, v = 10^-7
case = cases[1]
tapers[case] = [10, 100, 150, 250, 500, 1000, 2000]
extents[case] = [120, 120, 120, 240, 240, 240, 240]

# M = 5, v = 10^-6
case = cases[2]
tapers[case] = [2, 100, 200, 400]
extents[case] = [90, 120, 180, 240]

# M = 5, v = 10^-7
case = cases[3]
tapers[case] = [2, 200, 400, 800]
extents[case] = [120, 120, 180, 180]

# Other Cases #
case_offset = 1.08
taper_offset = 0.96

# M = 1.41, v = 10^-6
case = cases[4]
tapers[case] = np.array([100, 150])
extents[case] = [120, 240]

# M = 2, v = 10^-6
case = cases[5]
tapers[case] = np.array([100, 150, 250])
extents[case] = [120, 240, 240]

# M = 3.16, v = 10^-6
case = cases[6]
tapers[case] = np.array([100, 150, 250])
extents[case] = [120, 180, 240]

# M = 3.16, v = 4 * 10^-7
case = cases[7]
tapers[case] = np.array([100, 200])
extents[case] = [120, 180]

# M = 3.16, v = 1.33 * 10^-7
case = cases[8]
tapers[case] = np.array([200])
extents[case] = [120]

### For SVM later ###
concentrated_x = []
concentrated_y = []

elongated_x = []
elongated_y = []

### Read in model (assume it exists already)
f = open("model.dat", "r")
start = False

# Fit Parameters
w1_fit = 0
w2_fit = 0
b_fit = 0

for line in f:
	#print line[:-1]
	line_sp = line.split(" ")
	if start:
		alpha_y = float(line_sp[0])

		f1 = line_sp[1].split(":")
		x1 = float(f1[-1])

		f2 = line_sp[2].split(":")
		x2 = float(f2[-1])

		#print "x1, x2:", x1, x2

		# Update Model
		w1_fit += alpha_y * x1
		w2_fit += alpha_y * x2

	if "threshold" in line:
		b_fit = float(line_sp[0])
		start = True

print "w1: %.3f, w2: %.3f, b: %.3f" % (w1_fit, w2_fit, b_fit)

def fit(x):
	log_x = np.log(x) / np.log(10)

	m = -w1_fit / w2_fit
	b = b_fit / w2_fit

	print "m: %.3f, b: %.3f" % (m, 10**b)

	y =  np.power(10, m * log_x + b)
	print y
	#y = x**(m) * 10.0**(b)
	#print y

	return y


#### PLOTTING ####
fig = plot.figure()
ax = fig.add_subplot(1, 1, 1)

fontsize = 18

const = 90
power = 0.25
#plot.plot([0.5, 400], [const * 0.5**(power), const * 400**(power)], c = "g")
plot.plot([0.5, 400], [fit(0.5), fit(400)], c = "g")

for c, case in enumerate(cases):
	tapers_i = tapers[case]
	extents_i = extents[case]

	for i, (taper, extent) in enumerate(zip(tapers_i, extents_i)):
		if extent > 180:
			extent_type = "x"
			# For SVM
			elongated_x.append(np.log(case) / np.log(10))
			elongated_y.append(np.log(taper) / np.log(10))
		elif extent == 180:
			extent_type = "_"
		else:
			extent_type = "o"
			# For SVM
			concentrated_x.append(np.log(case) / np.log(10))
			concentrated_y.append(np.log(taper) / np.log(10))

		if c < 4:
			color = "blue"
			case_plot = 1.0 * case
			taper_plot = 1.0 * taper
		else:
			color = "red"
			case_plot = case_offset * case
			taper_plot = taper_offset * taper

		plot.scatter([case_plot], [taper_plot], marker = extent_type, c = color, s = 125)

#### Make Support Vector Machine ####

# Write Examples to File (in svm_light format)
f = open("train.dat", "w")
for (xi, yi) in zip(concentrated_x, concentrated_y):
	f.write("-1 1:%.3f 2:%.3f\n" % (xi, yi))

for (xi, yi) in zip(elongated_x, elongated_y):
	f.write("1 1:%.3f 2:%.3f\n" % (xi, yi))

f.close()

# Legend
legend_text = ""
legend_text += "        " + r"$\phi_{min}$ $>$ $180^{\circ}$" + "\n"
legend_text += "        " + r"$\phi_{min}$ $\approx$ $180^{\circ}$" + "\n"
legend_text += "        " + r"$\phi_{min}$ $<$ $180^{\circ}$"
plot.text(0.65, 1300, legend_text, fontsize = fontsize - 4, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1, pad = 7.0))

xtext = 0.78
ytext = 1355
legend_size = 90
color_sym = "darkblue"
plot.scatter(xtext, ytext * 10**(0.2), marker = "x", c = color_sym, s = legend_size)
plot.scatter(xtext, ytext * 10**(0.1), marker = "_", c = color_sym, s = legend_size)
plot.scatter(xtext, ytext, marker = "o", c = color_sym, s = legend_size)

# Annotate
plot.xlabel(r"$q^2 Re$ / $[3 \times 10^{-4}]$", fontsize = fontsize + 2)
plot.ylabel("Jupiter-Mass Growth Time (planet orbits)", fontsize = fontsize)
plot.title("Azimuthal Extents", fontsize = fontsize + 2)

# Axes
plot.xlim(0.5, 400)
plot.ylim(30, 3000)

plot.xscale("log")
plot.yscale("log")

axis_labels = ["50", "100", "200", "500", "1000", "2000"]
tick_locations = [int(x) for x in axis_labels]

ax.set_yticks(tick_locations)
ax.set_yticklabels(axis_labels)

# 2nd Axis
ax2 = ax.twinx()
ax2.set_yscale("log")
ax2.set_ylim(30**(-1), 3000**(-1))
#ax2.set_yticks()
#ax2.set_yticklabels(twin_axis_labels)
ax2.set_ylabel(r"$<\dot{q}>$ $\times$ $T_p$", rotation = 270, labelpad = 10, fontsize = fontsize + 2)

# Save, Show, and Close
plot.savefig("scattered_extents.png", bbox_inches = "tight")
plot.savefig("scattered_extents.pdf", bbox_inches = "tight", format = "pdf")
plot.show()
plot.cla()


