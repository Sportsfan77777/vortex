"""
compares variables for different taperings

Usage:
setupTriggerComparison.py
"""

import numpy as np
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc
from matplotlib import ticker as ticker

import pickle

### Data ###
# M1, v6
case1_x = 1.0 / np.array(pickle.load(open("case1_trigger_tapers.p", "rb"))) #[10, 250, 500, 1000]
case1_y = pickle.load(open("case1_triggers.p", "rb")) #[1000, 700, 400, 200]

# M1, v7
case2_x = 1.0 / np.array(pickle.load(open("case2_trigger_tapers.p", "rb"))) #[10, 500, 1000, 2000]
case2_y = pickle.load(open("case2_triggers.p", "rb")) #[7000, 1500, 1000, 600]

# M5, v6
case3_x = 5.0 / np.array(pickle.load(open("case3_trigger_tapers.p", "rb")))#[10, 500, 1000, 2000]
case3_y = pickle.load(open("case3_triggers.p", "rb")) #[1500, 1200, 500, 300]

# M5, v7
case4_x = 5.0 / np.array(pickle.load(open("case4_trigger_tapers.p", "rb")))#[10, 1000, 2000, 4000]
case4_y = pickle.load(open("case4_triggers.p", "rb")) #[14000, 10000, 8000, 5000]

#### Helper Function ####

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

#### PLOTTING ####
# Parameters
linewidth = 4
vertical_linewidth = 2
brace_linewidth = 2
fontsize = 17
labelsize = 15
markersize = 9
my_dpi = 100

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

max_value = 0.375
vertical_min_value = 10**(-5)
vertical_max_value = 3 * max_value

### Labels ###
case1_label = r"$1$ $M_J$, $\nu = 10^{-6}$"
case2_label = r"$1$ $M_J$, $\nu = 10^{-7}$"
case3_label = r"$5$ $M_J$, $\nu = 10^{-6}$"
case4_label = r"$5$ $M_J$, $\nu = 10^{-7}$"

### Setup Figure ###
figure = plot.figure(figsize = (850 / my_dpi, 600 / my_dpi), dpi = my_dpi)
ax = figure.add_subplot(1, 1, 1)

# Limits
log_x = True
log_y = False

if log_x:
    plot.xscale("log")
    plot.xlim(80**(-1), 2500**(-1))
    #ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
else:
    plot.xlim(-100, 2100)

if log_y:
    plot.yscale("log")
    max_value *= 2
    plot.ylim(0, max_value)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

    y_values = [1.0, 300, 1000, 3000, 10000]
    ax.set_yticks(y_values)
else:
    plot.ylim(0.125, max_value)

#x_values = [250, 500, 1000, 2000]
#ax.set_xticks(x_values)

### Vertical Lines ###
style4 = [5, 5]
line4_20 = 112 # M_dot = 10^-4 M_Jup / yr at 20 AU
line4_10 = 316 # M_dot = 10^-4 M_Jup / yr at 10 AU
line4_5 = 894 # M_dot = 10^-4 M_Jup / yr at 5 AU

style5 = [16, 10, 4, 10]
line5_20 = 1118 # M_dot = 10^-5 M_Jup / yr at 20 AU
line5_10 = 3162 # M_dot = 10^-5 M_Jup / yr at 10 AU

# Lines
#plot.plot([line4_20, line4_20], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, dashes = style4)
#plot.plot([line4_10, line4_10], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, dashes = style4)
#plot.plot([line4_5, line4_5], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, dashes = style4)#

#plot.plot([line5_20, line5_20], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, dashes = style5)
#plot.plot([line5_10, line5_10], [vertical_min_value, vertical_max_value], linewidth = vertical_linewidth, dashes = style5)

## Labels ##
# 10^-4
text_y4 = max_value + 0.005
#plot.text(line4_20, text_y4, "20", horizontalalignment = 'center')
#plot.text(line4_10, text_y4, "10", horizontalalignment = 'center')
#plot.text(line4_5, text_y4, "5", horizontalalignment = 'center')

midpoint4 = (line4_20 + line4_5) / 2.0
#plot.text(midpoint4, text_y4 + 0.03, r"$\dot{M}$ = $10^{-4}$ $M_J$ / yr" + "\nat 20, 10, and 5 AU", horizontalalignment = 'center')

# 10^-5
text_y5 = max_value + 0.010
#plot.text(line5_20, text_y5, "20", horizontalalignment = 'center')
#plot.text(line5_10, text_y5, "10", horizontalalignment = 'center')

midpoint5 = (line5_20 + line5_10) / 2.0
#plot.text(midpoint5, text_y5 + 0.03, r"$\dot{M}$ = $10^{-5}$ $M_J$ / yr at 20 and 10 AU", horizontalalignment = 'center')

# Braces
brace4_x, brace4_y = range_brace(line4_20, line4_5, height = 0.015)
#plot.plot(brace4_x, text_y4 + brace4_y + 0.01, color = "black", linewidth = brace_linewidth, clip_on = False)

brace5_x, brace5_y = range_brace(line5_20, line5_10, height = 0.015)
#plot.plot(brace5_x, text_y5 + brace5_y + 0.01, color = "black", linewidth = brace_linewidth, clip_on = False)

#### Curves ####
colors = [1, "orange", "gold", "deepskyblue", "darkblue"]
#colors = [1, "gold", "deepskyblue", "goldenrod", "darkblue"]

plot.plot(case3_x[1:], case3_y, marker = "p", c = colors[3], markersize = markersize + 1, linewidth = linewidth, label = case3_label) # M5, v6
plot.plot(case1_x[1:], case1_y, marker = "^", c = colors[1], markersize = markersize, linewidth = linewidth, label = case1_label) # M1, v6
plot.plot(case4_x[1:], case4_y, marker = "s", c = colors[4], markersize = markersize, linewidth = linewidth, label = case4_label) # M5, v7
plot.plot(case2_x[1:], case2_y, marker = "*", c = colors[2], markersize = markersize + 3, linewidth = linewidth, label = case2_label) # M1, v7

#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(round(x)))))

# Annotate
plot.xlabel(r"Averaged Accretion Rate ($M_J$ / planet orbit)", fontsize = fontsize)
plot.ylabel("Trigger Mass ($M_J$)", fontsize = fontsize)

legend = plot.legend(loc = "upper right") #, bbox_to_anchor = (legend_x, legend_y))

# Twin Axis
twin_axis_labels = ["100", "200", "500", "1000", "2000"]
twin_tick_locations = [1.0 / int(x) for x in twin_axis_labels]

print twin_tick_locations

ax2 = ax.twiny()
ax2.set_xscale("log")
ax2.set_xlim(80**(-1), 2500**(-1))
ax2.set_xticks(twin_tick_locations)
ax2.set_xticklabels(twin_axis_labels)
ax2.set_xlabel(r"Jupiter-Mass Growth Time (in planet orbits)", fontsize = fontsize, labelpad = 7)

# Save + Close
plot.savefig("trigger_comparison.png", bbox_inches = "tight")
plot.savefig("trigger_comparison.pdf", bbox_inches = "tight", format = "pdf")
plot.show()