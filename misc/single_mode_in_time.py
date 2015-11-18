"""
plots single mode
"""

import sys
import numpy as np 
from matplotlib import pyplot as plot

import subprocess

### Movie Commands ###
def make_movie(mode):
    # Movie Parameters
    fps = 5

    directory = "mode%d" % mode

    path = directory + "/" + directory + "_%03d.png"
    output = directory + "/" + directory + ".mov" 

    # Movie Command
    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    split_command = command.split()
    subprocess.Popen(split_command)


xs = np.linspace(0, 10, 100)
times = lambda max_t : np.linspace(0, max_t, 10 * max_t)

# Plotting
fontsize = 20
linewidth = 3

def make_plot(mode, t_max = 0):
    ts = times(t_max)
    for time in ts:
        """ sin wave with frequency 'mode' """
        ys = [np.sin(mode * x - np.pi * time / 10) for x in xs]

        figure = plot.figure(figsize = (10, 3))

        plot.plot(xs, ys, 'b', linewidth = linewidth + 1)

        # Limit
        plot.xlim(0, 2 * np.pi)
        plot.ylim(-1, 1)

        # Annotate
        plot.xlabel(r"$\phi$", fontsize = fontsize)
        plot.title("Mode $m = %d$" % mode, fontsize = fontsize + 2)

        # Save and Close
        plot.savefig("mode%d/mode%d_%03d.png" % (mode, mode, time), bbox_inches = 'tight')
        #plot.show()
        plot.close(figure)

mode_num = int(sys.argv[1])

make_plot(mode_num, t_max = 40)

make_movie(mode_num)


