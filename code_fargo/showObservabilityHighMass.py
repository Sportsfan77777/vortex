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
    return [arr[2], arr[3], arr[0], arr[1]]

### Data ###
# M5, v6
case3_x = pickle.load(open("case3_tapers.p", "rb")) #[10, 500, 1000, 2000]
case3_x = [case3_x[1], case3_x[-2]]

case3_y = pickle.load(open("case3_lifetimes.p", "rb")) #[1380, 1350, 660, 610]
case3_y = [case3_y[1], case3_y[-2]]

# M5, v7
case4_x = pickle.load(open("case4_tapers.p", "rb")) #[10, 1000, 2000, 4000]
case4_x = [case4_x[-2], case4_x[-1]]

case4_y = pickle.load(open("case4_lifetimes.p", "rb")) #[19990, 16730, 14680, 6770]
case4_y = [case4_y[-2], case4_y[-1]]

### Gather Cases ###
tapers = [case3_x, case4_x]
cases_y = [case3_y, case4_y]

### Convert Orbits to Years ###
xs = np.linspace(0.1, 100, 400)
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
dashes = [[10, 8], [4, 4], [10**5, 1]]
linewidth = 4

alpha_max = 0.35
soft_alphas = [alpha_max, alpha_max, alpha_max]

alpha_max = 0.9
hard_alphas = [alpha_max, alpha_max, alpha_max]

### Labels ###
case3_label = r"$\nu_{disk} = 10^{-6}$, $T_{growth} = "
case4_label = r"$\nu_{disk} = 10^{-7}$, $T_{growth} = "

labels = [case3_label, case4_label]

### Cutoffs ###
limit = 10**5

case3_cutoffs = [np.power(limit / case3_x[-2], 2.0 / 3.0), np.power(limit / case3_x[-1], 2.0 / 3.0)]
case4_cutoffs = [np.power(limit / case4_x[-2], 2.0 / 3.0), np.power(limit / case4_x[-1], 2.0 / 3.0)]

cutoffs = [case3_cutoffs, case4_cutoffs]

print cutoffs

high_limit = 5 * limit

case3_high_cutoffs = [np.power(high_limit / case3_x[-2], 2.0 / 3.0), np.power(high_limit / case3_x[-1], 2.0 / 3.0)]
case4_high_cutoffs = [np.power(high_limit / case4_x[-2], 2.0 / 3.0), np.power(high_limit / case4_x[-1], 2.0 / 3.0)]

high_cutoffs = [case3_high_cutoffs, case4_high_cutoffs]

print high_cutoffs

## Set up figure ##
fig = plot.figure()
ax = fig.add_subplot(111)

### Curves ###

# Definition #
plot.plot(xs, limit * np.ones(len(xs)), c = "r", linewidth = 2, zorder = 50)

# Data #

for i, case_y in enumerate(cases_yrs):
    for j, ys in enumerate(case_y):
        label = ""
        #if j == 0:
        #   label = labels[i] + str(tapers[i][j]) + "$"
        label = labels[i] + str(tapers[i][j]) + "$ $T_p$"
        print j, label

        cutoff = np.searchsorted(xs, cutoffs[i][j]) # convert cutoff from AU to index

        if tapers[i][j] == 10:
            plot.plot(xs, ys, c = colors[i], alpha = soft_alphas[j], linewidth = linewidth, dashes = dashes[j], label = label, zorder = 10)
        else:
            if len(case_y) == 2:
                j +=1 # Switch for dashes only

            plot.plot(xs[:cutoff], ys[:cutoff], c = colors[i], alpha = soft_alphas[j], linewidth = linewidth, dashes = dashes[j], zorder = 10)
            plot.plot(xs[cutoff:], ys[cutoff:], c = colors[i], alpha = hard_alphas[j], linewidth = linewidth + 1, dashes = dashes[j], label = label, zorder = 10)

            if len(case_y) == 2:
                j -= 1 # Switch for dashes only

        # Mark 10^5 years
        if tapers[i][j] > 10:
            cutoff_value = cutoffs[i][j]
            cutoff = np.searchsorted(xs, cutoff_value) # convert cutoff from AU to index
            plot.scatter(xs[cutoff], ys[cutoff], c = colors[i], marker = "s", s = 120, zorder = 100)

        # Mark 5 x 10^5 years
        if tapers[i][j] > 10:
            high_cutoff_value = high_cutoffs[i][j]
            high_cutoff = np.searchsorted(xs, high_cutoff_value) # convert cutoff from AU to index
            plot.scatter(xs[high_cutoff], ys[high_cutoff], c = colors[i], marker = "^", s = 120, zorder = 100)

# Axes
plot.xlim(0, 50)
plot.ylim(10**(2.5), 10**(6.7))
plot.yscale("log")

# Annotate
plot.xlabel("Planet Semimajor Axis (in AU)", fontsize = fontsize)
plot.ylabel("Vortex Lifetime (in years)", fontsize = fontsize)
plot.title(r"Vortex Observability: $m_p = 5$ $M_J$", fontsize = fontsize + 1, y = 1.02)

handles, labels = ax.get_legend_handles_labels()
plot.legend(re_sort(handles), re_sort(labels), loc = "lower right")

# Save + Show
plot.savefig("observability_highMass.png", bbox_inches = "tight")
plot.savefig("observability_highMass.pdf", bbox_inches = "tight", format = "pdf")
plot.show()



