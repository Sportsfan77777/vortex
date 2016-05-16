"""
compares variables for different taperings

Usage:
setupTaperComparison.py
"""

import numpy as np
from matplotlib import pyplot as plot
from matplotlib import ticker as ticker

import pickle

### Data ###
# M1, v6
case1_x = [10, 250, 500, 1000]
case1_y = [1000, 700, 400, 200]

# M1, v7
case2_x = [10, 500, 1000, 2000]
case2_y = [7000, 1500, 1000, 600]

# M5, v6
case3_x = [10, 500, 1000, 2000]
case3_y = [1500, 1200, 500, 300]

# M5, v7
case4_x = [10, 1000, 2000, 4000]
case4_y = [14000, 10000, 8000, 5000]

#### PLOTTING ####
# Parameters
linewidth = 4
vertical_linewidth = 2
fontsize = 14
markersize = 9

max_value = 15000
vertical_min_value = 1
vertical_max_value = 3 * max_value

### Labels ###
case1_label = r"$1 M_J$, $\nu = 10^{-6}$"
case2_label = r"$1 M_J$, $\nu = 10^{-7}$"
case3_label = r"$5 M_J$, $\nu = 10^{-6}$"
case4_label = r"$5 M_J$, $\nu = 10^{-7}$"

### Setup Figure ###
figure = plot.figure()
ax = figure.add_subplot(1, 1, 1)

### Vertical Lines ###
style4 = "--"
line4_20 = 112 # M_dot = 10^-4 M_Jup / yr at 20 AU
line4_10 = 316 # M_dot = 10^-4 M_Jup / yr at 10 AU
line4_5 = 894 # M_dot = 10^-4 M_Jup / yr at 5 AU

style5 = [16, 10, 4, 10]
line5_20 = 1118 # M_dot = 10^-5 M_Jup / yr at 20 AU
line5_10 = 3162 # M_dot = 10^-5 M_Jup / yr at 10 AU

plot.plot([line4_20, line4_20], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, linestyle = style4)
plot.plot([line4_10, line4_10], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, linestyle = style4)
plot.plot([line4_5, line4_5], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, linestyle = style4)

plot.plot([line5_20, line5_20], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, dashes = style5)
plot.plot([line5_10, line5_10], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, dashes = style5)

#### Curves ####
plot.plot(case4_x, case4_y, marker = "s", markersize = markersize, linewidth = linewidth, label = case4_label) # M5, v7
plot.plot(case2_x, case2_y, marker = "*", markersize = markersize + 3, linewidth = linewidth, label = case2_label) # M1, v7
plot.plot(case3_x, case3_y, marker = "p", markersize = markersize + 1, linewidth = linewidth, label = case3_label) # M5, v6
plot.plot(case1_x, case1_y, marker = "^", markersize = markersize, linewidth = linewidth, label = case1_label) # M1, v6

# Limits
log_x = False
log_y = True

if log_x:
    plot.xscale("log")
    plot.xlim(5, 5000)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
else:
    plot.xlim(-100, 4100)

if log_y:
    plot.yscale("log")
    plot.ylim(100, 2 * max_value)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

    y_values = [100, 300, 1000, 3000, 10000]
    ax.set_yticks(y_values)
else:
    plot.ylim(0, max_value)

x_values = [10, 500, 1000, 2000, 4000]
ax.set_xticks(x_values)


#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(round(x)))))

# Annotate
plot.xlabel("Taper Time (in planet orbits)", fontsize = fontsize)
plot.ylabel("Vortex Lifetime (in planet orbits)", fontsize = fontsize)

legend = plot.legend(loc = "lower right") #, bbox_to_anchor = (legend_x, legend_y))

# Save + Close
plot.savefig("lifetime_comparison.png")
plot.show()