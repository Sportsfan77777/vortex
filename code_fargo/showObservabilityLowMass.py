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
    #return [arr[2], arr[3], arr[4], arr[0], arr[1]]
    return arr

### Data ###
# M1, v6
case1_x = pickle.load(open("case1_tapers.p", "rb")) #[10, 250, 500]
case1_x = [case1_x[0], case1_x[-1]]

case1_y = pickle.load(open("case1_lifetimes.p", "rb")) #[530, 250, 230]
case1_y = [case1_y[0], case1_y[-1]]

# M1, v7
case2_x = pickle.load(open("case2_tapers.p", "rb")) #[10, 500, 1000, 2000]
case2_x = [case2_x[0], case2_x[-2], case2_x[-1]]

case2_y = pickle.load(open("case2_lifetimes.p", "rb")) #[7600, 2280, 1190, 850]
case2_y = [case2_y[0], case2_y[-2], case2_y[-1]]

### Gather Cases ###
tapers = [case1_x, case2_x]
cases_y = [case1_y, case2_y]

### Convert Orbits to Years ###
xs = np.linspace(0.1, 100, 800)
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

size = 100

colors = ["g", "b", "k", "y"]
linestyles = ["-", "-.", "--"]
dashes = [[10, 8], [4, 4], [10**5, 1]]
linewidth = 4

alpha_max = 0.35
soft_alphas = [alpha_max, alpha_max, alpha_max]

alpha_max = 0.9
hard_alphas = [alpha_max, alpha_max, alpha_max]

### Titles ###
titles = [r"$\nu = 10^{-6}$", r"$\nu = 10^{-7}$"]

### Labels ###
case1_label = r"$T_{growth} = "
case2_label = r"$T_{growth} = "

labels = [case1_label, case2_label]

### Cutoffs ###
limit = 10**5

case1_cutoffs = [np.power(limit / case1_x[-1], 2.0 / 3.0)]
case2_cutoffs = [np.power(limit / case2_x[-2], 2.0 / 3.0), np.power(limit / case2_x[-1], 2.0 / 3.0)]

cutoffs = [case1_cutoffs, case2_cutoffs]

print cutoffs

high_limit = 5 * limit

case1_high_cutoffs = [np.power(high_limit / case1_x[-1], 2.0 / 3.0)]
case2_high_cutoffs = [np.power(high_limit / case2_x[-2], 2.0 / 3.0), np.power(high_limit / case2_x[-1], 2.0 / 3.0)]

high_cutoffs = [case1_high_cutoffs, case2_high_cutoffs]

print high_cutoffs

## Set up figure ##
fig = plot.figure(figsize = (7, 8))
plot.subplots_adjust(wspace = 0.03, hspace = 0.03)

# Data #

for i, case_y in enumerate(cases_yrs):
    # One subplot for each
    ax = fig.add_subplot(2, 1, 2 - i)

    # Definition #
    plot.plot(xs, limit * np.ones(len(xs)), c = "r", linewidth = 2, zorder = 1)

    for j, ys in enumerate(case_y):
        label = ""
        #if j == 0:
        #   label = labels[i] + str(tapers[i][j]) + "$"
        print i, j
        print labels
        label = labels[i] + str(tapers[i][j]) + "$ $T_p$"
        print label

        cutoff = np.searchsorted(xs, cutoffs[i][j - 1]) # convert cutoff from AU to index

        if j == 0:
            plot.plot(xs, ys, c = colors[i], alpha = soft_alphas[j], linewidth = linewidth, dashes = dashes[j], label = label, zorder = 10)
        else:
            if len(case_y) == 2:
                k = 2 # Switch for dashes only
            else:
                k = j

            plot.plot(xs[:cutoff], ys[:cutoff], c = colors[i], alpha = soft_alphas[k], linewidth = linewidth, dashes = dashes[j], zorder = 10)
            plot.plot(xs[cutoff:], ys[cutoff:], c = colors[i], alpha = hard_alphas[k], linewidth = linewidth + 1, dashes = dashes[j], label = label, zorder = 10)

            if len(case_y) == 2:
                k = 1 # Switch for dashes only
            else:
                k = j

        # Mark 10^5 years
        if tapers[i][j] > 10:
            cutoff_value = cutoffs[i][j - 1]
            cutoff = np.searchsorted(xs, cutoff_value) # convert cutoff from AU to index
            plot.scatter(xs[cutoff], ys[cutoff], c = colors[i], marker = "s", s = size, zorder = 100)

        # Mark 5 x 10^5 years
        if tapers[i][j] > 10:
            high_cutoff_value = high_cutoffs[i][j - 1]
            high_cutoff = np.searchsorted(xs, high_cutoff_value) # convert cutoff from AU to index
            plot.scatter(xs[high_cutoff], ys[high_cutoff], c = colors[i], marker = "^", s = size, zorder = 100)

    # Axes
    plot.xlim(0, 50)
    plot.ylim(10**(2.5), 10**(6.7))
    plot.yscale("log")

    # Annotate
    if i == 0:
        plot.xlabel("Planet Semimajor Axis (in AU)", fontsize = fontsize)
    else:
        ax.set_xticklabels([])
    plot.ylabel("Vortex Lifetime (in years)", fontsize = fontsize)
    plot.title(titles[i], fontsize = fontsize + 1, x = 0.11, y = 0.85, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))

    #handles, labels = ax.get_legend_handles_labels()
    #plot.legend(re_sort(handles), re_sort(labels), loc = "lower right")
    plot.legend(loc = "lower right")

    ## Setup Accretion Rate Legend ##
    xtext = 10
    margin = 1
    offset = 0

    accretion_text = r"      $\dot{M} = 10^{-5} M_J \rm{ / yr}$      "
    if i == 1:
        accretion_text += "\n" + r"      $\dot{M} = 2 \times 10^{-6} M_J \rm{ / yr}$"
        offset = 0.4

    plot.text(xtext, 10**(2.8), accretion_text, fontsize = fontsize - 4, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1, pad = 7.0))
    plot.scatter(xtext + margin, 10**(2.9 + offset), c = colors[i], marker = "s", s = size)
    if i == 1:
        plot.scatter(xtext + margin, 10**(2.9), c = colors[i], marker = "^", s = size)

plot.suptitle(r"Vortex Observability: $m_p = 1$ $M_J$", fontsize = fontsize + 1, y = 0.95)

# Save + Show
plot.savefig("observability_lowMass.png", bbox_inches = "tight")
plot.savefig("observability_lowMass.pdf", bbox_inches = "tight", format = "pdf")
plot.show()



