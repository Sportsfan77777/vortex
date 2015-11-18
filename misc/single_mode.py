"""
plots single mode
"""

import sys
import numpy as np 
from matplotlib import pyplot as plot
from matplotlib import ticker

import subprocess

xs = np.linspace(0, 2*np.pi, 100)
xs_deg = np.linspace(0, 360, 100)

ticks = ticker.MultipleLocator(base = 60)

# Plotting
fontsize = 20
linewidth = 3

def make_plot(mode):
    """ sin wave with frequency 'mode' """
    ys = [np.sin(mode * x) for x in xs]

    figure, ax = plot.subplots(figsize = (10, 3))

    plot.plot(xs_deg, ys, 'b', linewidth = linewidth + 1)

    # Limit
    plot.xlim(0, 360)
    plot.ylim(-1, 1)

    # Annotate
    plot.xlabel(r"$\phi$ (degrees)", fontsize = fontsize)
    ax.xaxis.set_major_locator(ticks)
    plot.title("Mode $m = %d$" % mode, fontsize = fontsize + 2)

    # Save and Close
    plot.savefig("mode%d.png" % (mode), bbox_inches = 'tight')
    #plot.show()
    plot.close(figure)

mode_num = int(sys.argv[1])

make_plot(mode_num)
