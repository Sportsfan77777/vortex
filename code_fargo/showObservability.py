"""
marks when vortex lifetime corresponds to above a threshold number of years

Usage:
showObservability.py
"""

import numpy as np
from matplotlib import pyplot as plot
from matplotlib import ticker as ticker
from matplotlib import rcParams as rc

import pickle

## Helper Function ##
def re_sort(arr):
	return [arr[3], arr[1], arr[2], arr[0]]

### Data ###
# M1, v6
case1_x = pickle.load(open("case1_tapers.p", "rb")) #[10, 250, 500]
#case1_x = [case1_x[0], case1_x[-1]]

case1_y = pickle.load(open("case1_lifetimes.p", "rb")) #[530, 250, 230]
#case1_y = [case1_y[0], case1_y[-1]]

# M1, v7
case2_x = pickle.load(open("case2_tapers.p", "rb")) #[10, 500, 1000, 2000]
case2_x = [case2_x[0], case2_x[-2], case2_x[-1]]

case2_y = pickle.load(open("case2_lifetimes.p", "rb")) #[7600, 2280, 1190, 850]
case2_y = [case2_y[0], case2_y[-2], case2_y[-1]]

# M5, v6
case3_x = pickle.load(open("case3_tapers.p", "rb")) #[10, 500, 1000, 2000]
case3_x = [case3_x[0], case3_x[-2], case3_x[-2]]

case3_y = pickle.load(open("case3_lifetimes.p", "rb")) #[1380, 1350, 660, 610]
case3_y = [case3_y[0], case3_y[-2], case3_y[-2]]

# M5, v7
case4_x = pickle.load(open("case4_tapers.p", "rb")) #[10, 1000, 2000, 4000]
case4_x = [case4_x[0], case4_x[-2], case4_x[-1]]

case4_y = pickle.load(open("case4_lifetimes.p", "rb")) #[19990, 16730, 14680, 6770]
case4_y = [case4_y[0], case4_y[-2], case4_y[-1]]

### Gather Cases ###
tapers = [case1_x, case2_x, case3_x, case4_x]
cases_y = [case1_y, case2_y, case3_y, case4_y]

### Convert Orbits to Years ###
xs = np.linspace(0.1, 50, 400)
Ts = np.array([x**(1.5) for x in xs])

cases_yrs = []
for case_y in cases_y:
	cases_yrs.append([Ts * lifetime for lifetime in case_y])

### PLOTTING ###
# Parameters #
fontsize = 16

labelsize = 14
rc["xtick.labelsize"] = labelsize
rc["ytick.labelsize"] = labelsize

colors = ["g", "b", "k", "y"]
linestyles = ["-", "-.", "--"]
dashes = [[12, 10], [3, 1], [10**5, 1]]
linewidth = 4

alpha_min = 0.15
alpha_max = 0.35
soft_alphas = [alpha_max, alpha_min, alpha_max]

alpha_min = 0.35
alpha_max = 0.9
hard_alphas = [alpha_max, alpha_min, alpha_max]

### Labels ###
case1_label = r"$1$ $M_J$, $\nu = 10^{-6}$, $T_{growth} = "
case2_label = r"$1$ $M_J$, $\nu = 10^{-7}$, $T_{growth} = "
case3_label = r"$5$ $M_J$, $\nu = 10^{-6}$, $T_{growth} = "
case4_label = r"$5$ $M_J$, $\nu = 10^{-7}$, $T_{growth} = "

labels = [case1_label, case2_label, case3_label, case4_label]

### Cutoffs ###
limit = 10**5

case1_cutoff = np.power(limit / case1_x[-1], 2.0 / 3.0)
case2_cutoff = np.power(limit / case2_x[-1], 2.0 / 3.0)
case3_cutoff = np.power(limit / case3_x[-1], 2.0 / 3.0)
case4_cutoff = np.power(limit / case4_x[-1], 2.0 / 3.0)

cutoffs = [case1_cutoff, case2_cutoff, case3_cutoff, case4_cutoff]

print cutoffs

## Set up figure ##
fig = plot.figure()
ax = fig.add_subplot(111)

### Curves ###

# Definition #
plot.plot(xs, limit * np.ones(len(xs)), c = "r", linewidth = 2)

# Data #

for i, case_y in enumerate(cases_yrs):
	for j, ys in enumerate(case_y):
		label = ""
		#if j == 0:
		#	label = labels[i] + str(tapers[i][j]) + "$"
		if j == len(case_y) - 1:
			label = labels[i] + str(tapers[i][j]) + "$"

		cutoff = np.searchsorted(xs, cutoffs[i]) # convert cutoff from AU to index

		if j == 0:
			plot.plot(xs, ys, c = colors[i], alpha = soft_alphas[j], linewidth = linewidth, dashes = dashes[j], label = label)
		if j == 2:
			plot.plot(xs[:cutoff], ys[:cutoff], c = colors[i], alpha = soft_alphas[j], linewidth = linewidth, dashes = dashes[j])
			plot.plot(xs[cutoff:], ys[cutoff:], c = colors[i], alpha = hard_alphas[j], linewidth = linewidth + 1, dashes = dashes[j], label = label)

# Axes
plot.xlim(0, xs[-1])
plot.ylim(10**(2.5), 10**(6.5))
plot.yscale("log")

# Annotate
plot.xlabel("Planet Semimajor Axis (in AU)", fontsize = fontsize)
plot.ylabel("Vortex Lifetime (in years)", fontsize = fontsize)
plot.title("Vortex Observability", fontsize = fontsize + 1)

handles, labels = ax.get_legend_handles_labels()
plot.legend(re_sort(handles), re_sort(labels), loc = "lower right")

# Save + Show
plot.savefig("observability.png", bbox_inches = "tight")
plot.savefig("observability.pdf", bbox_inches = "tight", format = "pdf")
plot.show()



