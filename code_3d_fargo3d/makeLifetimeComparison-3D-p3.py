import numpy as np
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc
from matplotlib import ticker as ticker

from itertools import groupby
import pickle

###############################################################################

master_directories = []
lifetimes = []
final_masses = []

#master_directories.append(["h08_nu7_a167", "h08_nu7_a05", "h08_nu7_a02", "h08_nu7_a01"])
#lifetimes.append([3573, 3875, 7249, 9625])
#final_masses.append([1.22, 0.59, 0.33, 0.18])

#### Note: Re-do ALL of these!!! ####
### 2-D (2021) ###
master_directories.append(["h06_nu7_a50", "h06_nu7_a167", "h06_nu7_a05", "h06_nu7_a02"])
lifetimes.append([6712, 5822, 2236, 2424, 3900])
final_masses.append([0.21, 0.36, 0.63, 0.92, 1.05])

### 2-D low-res ###
master_directories.append(["h06_nu7_a01", "h06_nu7_a02", "h06_nu7_a04", "h06_nu7_a10", "h06_nu7_a167", "h06_nu7_a30", "h06_nu7_a70"])
lifetimes.append([10000, 10000, 7350, 6310, 6370, 9440, 10000])
final_masses.append([0.17, 0.26, 0.37, 0.61, 0.80, 1.04, 1.35])

### 3-D low-res ###
master_directories.append(["h06_nu7_a01", "h06_nu7_a02", "h06_nu7_a04", "h06_nu7_a10", "h06_nu7_a167", "h06_nu7_a30", "h06_nu7_a70"])
lifetimes.append([4505, 3990, 3335, 2690, 3045, 4000, 4795])
final_masses.append([0.22, 0.32, 0.44, 0.69, 0.89, 1.18, 1.51])


#master_directories.append(["h04_nu7_a100", "h04_nu7_a50", "h04_nu7_a167", "h04_nu7_a05"])
#lifetimes.append([2038, 1155, 1794, 2695])
#final_masses.append([0.67, 0.53, 0.36, 0.21])

#master_directories.append(["h08_nu6_a167", "h08_nu6_a05", "h08_nu6_a02"])
#lifetimes.append([1440, 1526, 0])
#final_masses.append([1.26, 0.70, 0.37])

#master_directories.append(["h06_nu6_a50", "h06_nu6_a167", "h06_nu6_a05"])
#lifetimes.append([559, 1089, 930])
#final_masses.append([1.18, 0.84, 0.52])

###############################################################################

### Helper Functions ###

def get_final_mass(directory):
    planet_data = np.loadtxt("../%s/planet0.dat" % directory)
    times = planet_data[:, 0]
    masses = planet_data[:, 7] + planet_data[:, 8] # Core Mass + Accreted Mass

    time_index = np.searchsorted(times, 2000)
    final_mass = masses[time_index] / 0.001
    return final_mass


### PLOTTING ###
linewidth = 6
labelsize = 17
fontsize = 20
markersize = 16
dpi = 100
alpha = 0.8

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

linestyles = ["-", "-", "-", "-", "--"]
linewidths = [linewidth, linewidth, linewidth, linewidth]
#colors = ["darkblue", "gold", "deepskyblue", "orange", "k"]
#colors = ["b", "#17becf", "gold", "orange", "r"]
colors = ["b", "orange", "gold"]
#markers = ["s", "*", "p", "^", "h"]
markers = ["s", "^", "*", "h"]
labels = ["2-D (MH+ 2021)", "2-D low res", "3-D low-res"]

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (8, 7), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data
    for i, (directories, lifetimes_i, final_masses_i) in enumerate(zip(master_directories, lifetimes, final_masses)):
        scale_height = float(directories[0].split("_")[0][1:]) / 100.0
        log_viscosity = float(directories[0].split("_")[1][2:]) - 2.0

        alpha_coefficent = "3"
        if scale_height == 0.08:
            alpha_coefficent = "1.5"
        elif scale_height == 0.04:
            alpha_coefficent = "6"
        #label = r"$h =$ $%.02f$, $\alpha = %s \times 10^{-%d}$" % (scale_height, alpha_coefficent, log_viscosity)
        #label = r"$h =$ $%.02f$" % (scale_height)
        label = labels[i]

        xs = final_masses_i
        ys = lifetimes_i

        plot.plot(xs, ys, c = colors[i], marker = markers[i], markersize = markersize, linewidth = linewidths[i], linestyle = linestyles[i], label = label, zorder = 99, alpha = alpha)

        # Indicate vortices that aren't dead yet
        dx = 0; dy = 600; head_width = 0.05; head_length = 200

        for (xi, yi) in zip(xs, ys):
            if yi >= 10000 or yi == 3900:
                plot.arrow(xi, yi, dx, dy, linewidth = linewidths[i] - 1, head_width = head_width, head_length = head_length, fc = colors[i],  ec = colors[i], alpha = alpha)

        offset_x = 0.025; offset_y = 200
        #plot.text(xs[0] + offset_x, ys[0] + offset_y, "S", color = colors[i], fontsize = fontsize - 3)
        #plot.text(xs[1] + offset_x, ys[1] + offset_y, "S", color = colors[i], fontsize = fontsize - 3)

        # Create break in legend
        #if directories == master_directories[2]:
        #    plot.plot([1, 1], [1, 1], c = "white", label = " ")

    plot.legend(loc = 'lower right', fontsize = fontsize - 3)

    # Axes
    plot.xlim(0, 1.6)
    plot.ylim(0, 11000)

    # Annotate
    plot.xlabel(r"Planet Mass [Jupiter masses]", fontsize = fontsize)
    plot.ylabel(r"Lifetime [planet orbits]", fontsize = fontsize)
    plot.title("Total Lifetimes of Dust Asymmetry\n(w/ varied planet accretion rates) ", y = 1.02, fontsize = fontsize + 1)

    xtext = 0.25 * (plot.xlim()[1] - plot.xlim()[0]) + plot.xlim()[0]
    ytext = 0.90 * (plot.ylim()[1] - plot.ylim()[0]) + plot.ylim()[0]
    enter = 0.06 * (plot.ylim()[1] - plot.ylim()[0]) + plot.ylim()[0]

    text1 = r"$h = 0.06$"; text2 = r"$\alpha = 3 \times 10^{-5}$"
    plot.text(xtext, ytext, text1, fontsize = fontsize)
    plot.text(xtext, ytext - enter, text2, fontsize = fontsize)

    # Save + Show
    plot.savefig("present_lifetime_comparison-3D-p3.png", bbox_inches = "tight")
    plot.savefig("present_lifetime_comparison-3D-p3.pdf", bbox_inches = "tight", format = "pdf")
    plot.show()

make_plot()
