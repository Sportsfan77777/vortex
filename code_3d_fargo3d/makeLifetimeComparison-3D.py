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
### 2-D ###
master_directories.append(["h06_nu7_s0347", "h06_nu7_s0578", "h06_nu7_s0694", "h06_nu7_s0926", "h06_nu7_s1157", "h06_nu7_s1736", "h06_nu7_s2315", "h06_nu7_s3472"])
lifetimes.append([0, 900, 1170, 1170, 2920, 7500, 7500, 7500, 7500])
final_masses.append([0.10, 0.15, 0.17, 0.21, 0.29, 0.42, 0.51, 0.65, 0.91])

### 2-D migrating ###
master_directories.append(["h06_nu7_s1157", "h06_nu7_s1736", "h06_nu7_s2315", "h06_nu7_s3472", "h06_nu7_s5787"])
lifetimes.append([1760, 2620, 4120, 4875, 7500])
final_masses.append([0.29, 0.42, 0.51, 0.65, 0.91])

### 2-D low-res ###
master_directories.append(["h06_nu7_s0578", "h06_nu7_s0694", "h06_nu7_s0926", "h06_nu7_s1157", "h06_nu7_s1736", "h06_nu7_s2315", "h06_nu7_s3472", "h06_nu7_s4629", "h06_nu7_s5787"])
lifetimes.append([0, 1360, 1970, 2370, 2950, 3340, 3780, 5690, 7500])
final_masses.append([0.15, 0.17, 0.23, 0.28, 0.39, 0.50, 0.71, 0.94, 1.19]) # 1.19

### 3-D low-res ###
master_directories.append(["h06_nu7_s0694", "h06_nu7_s0809", "h06_nu7_s0926", "h06_nu7_s1157", "h06_nu7_s1736", "h06_nu7_s2315", "h06_nu7_s3472", "h06_nu7_s5787"])
lifetimes.append([0, 1300, 1830, 2710, 3140, 2990, 3290, 6000, 6000])
final_masses.append([0.15, 0.17, 0.20, 0.26, 0.40, 0.51, 0.80, 1.13, 1.47]) # 1.47


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

linestyles = ["-", "--", "-", "-", "--"]
linewidths = [linewidth, linewidth, linewidth, linewidth]
#colors = ["darkblue", "gold", "deepskyblue", "orange", "k"]
#colors = ["b", "#17becf", "gold", "orange", "r"]
colors = ["b", "b", "orange", "gold"]
#markers = ["s", "*", "p", "^", "h"]
markers = ["s", "p", "^", "*", "h"]
labels = ["2-D", "2-D w/ dust cutoff", "2-D low res", "3-D low-res"]

def make_plot():
    # Set up figure
    fig = plot.figure(figsize = (8, 6), dpi = dpi)
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
            if yi >= 6000:
                plot.arrow(xi, yi, dx, dy, linewidth = linewidths[i] - 1, head_width = head_width, head_length = head_length, fc = colors[i],  ec = colors[i], alpha = alpha)

        offset_x = 0.025; offset_y = 200
        #plot.text(xs[0] + offset_x, ys[0] + offset_y, "S", color = colors[i], fontsize = fontsize - 3)
        #plot.text(xs[1] + offset_x, ys[1] + offset_y, "S", color = colors[i], fontsize = fontsize - 3)

        # Create break in legend
        #if directories == master_directories[2]:
        #    plot.plot([1, 1], [1, 1], c = "white", label = " ")

    plot.legend(loc = 'lower right', fontsize = fontsize - 3)

    # Axes
    plot.xlim(0, 1.5)
    plot.ylim(0, 9000)

    # Annotate
    plot.xlabel(r"Planet Mass [Jupiter masses]", fontsize = fontsize)
    plot.ylabel(r"Lifetime [planet orbits]", fontsize = fontsize)
    plot.title("Total Lifetimes of Dust Asymmetry", y = 1.02, fontsize = fontsize + 1)

    xtext = 0.025 * (plot.xlim()[1] - plot.xlim()[0]) + plot.xlim()[0]
    ytext = 0.90 * (plot.ylim()[1] - plot.ylim()[0]) + plot.ylim()[0]
    enter = 0.08 * (plot.ylim()[1] - plot.ylim()[0]) + plot.ylim()[0]

    text1 = r"$h = 0.06$"; text2 = r"$\alpha = 3 \times 10^{-5}$"
    plot.text(xtext, ytext, text1, fontsize = fontsize)
    plot.text(xtext, ytext - enter, text2, fontsize = fontsize)

    # Save + Show
    plot.savefig("present_lifetime_comparison-3D.png", bbox_inches = "tight")
    plot.savefig("present_lifetime_comparison-3D.pdf", bbox_inches = "tight", format = "pdf")
    plot.show()

make_plot()
