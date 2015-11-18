"""
plots single mode
"""

import sys
import numpy as np 
from matplotlib import pyplot as plot
from matplotlib import ticker

import subprocess

### Movie Commands ###
def make_movie(mode):
    # Movie Parameters
    fps = 5

    directory = "growing_mode%d" % mode

    path = directory + "/" + directory + "_%03d.png"
    output = directory + "/" + directory + ".mov" 

    # Movie Command
    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    split_command = command.split()
    subprocess.Popen(split_command)


xs = np.linspace(0, 2 * np.pi, 100)
xs_deg = np.linspace(0, 360, 100)
times = lambda max_t : np.linspace(0, max_t, 10 * max_t)


def exp(r, x):
    return np.exp(r * x)

rate = 0.05
grow = lambda x : exp(rate, x)

# Plotting
fontsize = 20
linewidth = 3

ticks = ticker.MultipleLocator(base = 60)

def make_plot(mode, t_max = 0):
    ts = times(t_max + 1)
    for time in ts:
        """ sin wave with frequency 'mode' """
        ys = [np.sin(mode * x - np.pi * time / 10) * grow(time) for x in xs]

        figure, ax = plot.subplots(figsize = (10, 3))

        plot.plot(xs_deg, ys, 'b', linewidth = linewidth + 1)

        # Limit
        plot.xlim(0, 360)
        plot.ylim(-8, 8)

        # Annotate
        plot.xlabel(r"$\phi$ (degrees)", fontsize = fontsize)
        ax.xaxis.set_major_locator(ticks)
        plot.title("Mode $m = %d$" % mode, fontsize = fontsize + 2)

        # Save and Close
        plot.savefig("growing_mode%d/growing_mode%d_%03d.png" % (mode, mode, time), bbox_inches = 'tight')
        #plot.show()
        plot.close(figure)

mode_num = int(sys.argv[1])

make_plot(mode_num, t_max = 40)

make_movie(mode_num)


