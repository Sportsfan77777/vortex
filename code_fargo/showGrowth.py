"""
compares variables for different taperings

Usage:
showGrowth.py
"""

import numpy as np
from matplotlib import pyplot as plot
from matplotlib import ticker as ticker
from matplotlib import rcParams as rc

import pickle

growth_times = np.array([250, 500, 1000, 2000, 4000])

### Convert Orbits to Years ###
xs = np.linspace(0.1, 50, 400)
Ts = np.array([x**(1.5) for x in xs])

years = []
for time in growth_times:
	years.append(Ts * time)


### PLOTTING ###
# Parameters
fontsize = 16
linewidth = 4

labelsize = 14
rc["xtick.labelsize"] = labelsize
rc["ytick.labelsize"] = labelsize

colors = ["orange", "g", "k", "b", "y"]

## Set up figure ##
fig = plot.figure()
ax = fig.add_subplot(111)

# Definition #
limit = 10**5
plot.plot(xs, limit * np.ones(len(xs)), c = "r", linewidth = 2)
plot.plot(xs, 5 * limit * np.ones(len(xs)), c = "r", linewidth = 2)

# Curves
for i, (time, ys) in enumerate(zip(growth_times, years)):
	plot.plot(xs, ys, linewidth = linewidth, c = colors[i], label = r"$T_{growth} = " + str(time) + "$ $T_{planet}$")

# Axes
plot.xlim(0, xs[-1])
plot.ylim(10**(2.5), 10**(6.5))
plot.yscale("log")

# Annotate
plot.xlabel("Planet Semimajor Axis (in AU)", fontsize = fontsize)
plot.ylabel(r"$T_{growth}$ (in years)", fontsize = fontsize)
plot.title("", fontsize = fontsize + 1)

handles, labels = ax.get_legend_handles_labels()
plot.legend(reversed(handles), reversed(labels), loc = "lower right")

# Save + Show
plot.savefig("growth_in_years.png", bbox_inches = "tight")
plot.savefig("growth_in_years.pdf", bbox_inches = "tight", format = "pdf")
plot.show()

