"""
Create stacked bar chart counting peaks at different resolutions
"""

import numpy as np
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc

### Parameters ###

scale_height = 0.06
surface_density_zero = 3.472e-4

### Data ###

# Peak distributions (1, 2, 3, 4+) from t = 100 to 1600 (every 10)
beam_sizes = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0])
labels = ["1", "2", "3", "4+"]

peak_distributions = []
#peak_distributions.append([0.25, 0.33, 0.19, 0.23]) # 4
peak_distributions.append([37, 53, 26, 35]) # 5
#peak_distributions.append([0.25, 0.34, 0.21, 0.20]) # 6
peak_distributions.append([39, 55, 45, 12]) # 10
peak_distributions.append([45, 74, 30, 2]) # 15
peak_distributions.append([63, 51, 30, 7]) # 20
peak_distributions.append([64, 46, 26, 15]) # 25
peak_distributions.append([77, 54, 16, 4]) # 30

peak_distributions = np.array(peak_distributions)
num_frames = 151.0

### PLOTTING ###
linewidth = 6
labelsize = 17
dpi = 100

fontsize = 18
markersize = 14

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

linestyles = ["-", "-", "-", "--", "--"]
#colors = ["darkblue", "gold", "deepskyblue", "orange", "k"]
colors = ["b", "#17becf", "gold", "orange", "r"]
#markers = ["s", "*", "p", "^", "h"]
markers = ["s", "p", "*", "^", "h"]

def make_plot():
    fig = plot.figure(figsize = (7, 5), dpi = dpi)
    ax = fig.add_subplot(111)

    bottom = np.zeros(6)

    for i in range(len(peak_distributions[0, :])):
        peak_distribution = peak_distributions[:, i] / num_frames
        if i > 0:
            bottom += peak_distributions[:, i - 1] / num_frames
        plot.barh(beam_sizes, peak_distribution, color = colors[i], left = bottom, height = 3, label = labels[i])

    legend = plot.legend(loc = 'upper left', fontsize = fontsize - 4, facecolor = 'white', framealpha = 1.0)
    legend.set_zorder(150)

    # Axes
    plot.xlim(0, 1)

    # Annotate
    plot.xlabel(r"Peak Count Distribution", fontsize = fontsize)
    plot.ylabel(r"Beam Diameters [$^{\prime\prime}$]", fontsize = fontsize)

    ### Add times used and legend

    title = r"$t = 100$ to $1600$ $\mathrm{orbits}}$ [$m_\mathrm{p}(t)\ =\ 0.57$ $M_\mathrm{Jup}$]"
    plot.title("%s" % (title), y = 1.02, fontsize = fontsize + 3)

    title2 = r'$h = %.2f$   $\Sigma = %.3e$  (2-D)' % (scale_height, surface_density_zero)
    x_mid = 0.5; y_text = 1.18
    plot.text(x_mid, y_text * plot.ylim()[-1], title2, fontsize = fontsize + 3, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))

    # Save + Show
    scale_height_name = 100 * scale_height
    surface_density_name = int(round(1e7 * surface_density_zero, 0))
    plot.savefig("peak_counts_across_resolution-h%02d-s%04d.png" % (scale_height_name, surface_density_name), bbox_inches = "tight")
    plot.savefig("peak_counts_across_resolution-h%02d-s%04d.pdf" % (scale_height_name, surface_density_name), bbox_inches = "tight", format = "pdf")
    plot.show()


make_plot()