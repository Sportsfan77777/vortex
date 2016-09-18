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
tapers[case] = np.array([100, 200, 400])
extents[case] = [120, 180, 180]

# M = 3.16, v = 1.33 * 10^-7
case = cases[8]
tapers[case] = np.array([200])
extents[case] = [120]

### For SVM later ###
concentrated_x = []
concentrated_y = []

elongated_x = []
elongated_y = []

## Helper Functions ##

def convert_to_log(arr):
	range_arr = (arr - arr[0]) / (arr[0] + arr[-1])
	print range_arr
	new_range = arr[0] * np.exp(range_arr * (np.log(arr[0]) + np.log(arr[-1])))
	print new_range
	return new_range

def range_brace(x_min, x_max, mid=0.75, 
                beta1=50.0, beta2=100.0, height=1, 
                initial_divisions=11, resolution_factor=1.5):
    # Source: http://stackoverflow.com/questions/1289681/drawing-braces-with-pyx
    # determine x0 adaptively values using second derivitive
    # could be replaced with less snazzy:
    #   x0 = np.arange(0, 0.5, .001)
    x0 = np.array(())
    tmpx = np.linspace(0, 0.5, initial_divisions)
    tmp = beta1**2 * (np.exp(beta1*tmpx)) * (1-np.exp(beta1*tmpx)) / np.power((1+np.exp(beta1*tmpx)),3)
    tmp += beta2**2 * (np.exp(beta2*(tmpx-0.5))) * (1-np.exp(beta2*(tmpx-0.5))) / np.power((1+np.exp(beta2*(tmpx-0.5))),3)
    for i in range(0, len(tmpx)-1):
        t = int(np.ceil(resolution_factor*max(np.abs(tmp[i:i+2]))/float(initial_divisions)))
        x0 = np.append(x0, np.linspace(tmpx[i],tmpx[i+1],t))
    x0 = np.sort(np.unique(x0)) # sort and remove dups
    # half brace using sum of two logistic functions
    y0 = mid*2*((1/(1.+np.exp(-1*beta1*x0)))-0.5)
    y0 += (1-mid)*2*(1/(1.+np.exp(-1*beta2*(x0-0.5))))
    # concat and scale x
    x = np.concatenate((x0, 1-x0[::-1])) * float((x_max-x_min)) + x_min
    y = np.concatenate((y0, y0[::-1])) * float(height)
    return (x,y)

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
legend_text += "        " + r"$\left[\Delta \phi\right]_\mathrm{min}$ $>$ $180^{\circ}$" + "\n"
legend_text += "        " + r"$\left[\Delta \phi\right]_\mathrm{min}$ $\approx$ $180^{\circ}$" + "\n"
legend_text += "        " + r"$\left[\Delta \phi\right]_\mathrm{min}$ $<$ $180^{\circ}$"
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

### Add Braces to Mark Red Points ###
top_brace_red = 70
bottom_brace_blue = 580

red_text_y = 38
blue_text_y = 850

# (1)

start = case_offset * cases[4]
end = case_offset * cases[6] * 0.97

brace1_x, brace1_y = range_brace(start, end, height = 20)
log_brace1_x = np.logspace(np.log10(start), np.log10(end), len(brace1_x))
midpoint = np.median(log_brace1_x)

plot.plot(log_brace1_x, top_brace_red - brace1_y, color = "r", linewidth = 2)
plot.text(midpoint, red_text_y, r"$\nu = 10^{-6}$", color = "r", horizontalalignment = 'center', fontsize = fontsize - 2)

# (2)

start = case_offset * cases[-3] * 1.03
end = case_offset * cases[-1]

brace2_x, brace2_y = range_brace(start, end, height = 20)
log_brace2_x = np.logspace(np.log10(start), np.log10(end), len(brace2_x))
midpoint = np.median(log_brace2_x)

plot.plot(log_brace2_x, top_brace_red - brace2_y, color = "r", linewidth = 2)
plot.text(midpoint, red_text_y, r"$m_p = 3.16$ $M_J$", color = "r", horizontalalignment = 'center', fontsize = fontsize - 2)

# (3)

start = case_offset * cases[0] * 1.1
end = case_offset * cases[1] * 0.8

brace3_x, brace3_y = range_brace(start, end, height = 150)
log_brace3_x = np.logspace(np.log10(start), np.log10(end), len(brace2_x))
midpoint = np.median(log_brace3_x)

plot.plot(log_brace3_x, bottom_brace_blue + brace3_y, color = "b", linewidth = 2)
plot.text(midpoint, blue_text_y, r"$m_p = 1$ $M_J$", color = "b", horizontalalignment = 'center', fontsize = fontsize - 2)

# (4)

start = case_offset * cases[2]
end = case_offset * cases[3] * 0.85

brace4_x, brace4_y = range_brace(start, end, height = 400)
log_brace4_x = np.logspace(np.log10(start), np.log10(end), len(brace2_x))
midpoint = np.median(log_brace4_x)

plot.plot(log_brace4_x, bottom_brace_blue + 500 + brace4_y, color = "b", linewidth = 2)
plot.text(midpoint, blue_text_y + 800, r"$m_p = 5$ $M_J$", color = "b", horizontalalignment = 'center', fontsize = fontsize - 2)


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


