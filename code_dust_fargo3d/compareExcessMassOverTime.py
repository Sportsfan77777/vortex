"""
compares excess masses with different parameters

Usage:
python compareExcessMassesOverTime.py
"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot
from matplotlib import gridspec

from itertools import groupby

from pylab import rcParams # replace with rc ???
from pylab import fromfile

import util
from readTitle import readTitle

###############################################################################

#directories = ["h08_nu7_a167", "h08_nu7_a05", "h08_nu7_a02", "h08_nu7_a01"]
#directories = ["h06_nu6_a50", "h06_nu6_a167", "h06_nu6_a05", "h06_nu6_a02"]
directories = ["h06_nu7_a50", "h06_nu7_a167", "h06_nu7_a05", "h06_nu7_a02"]
#directories = ["h06_nu0_a50", "h06_nu0_a167", "h06_nu0_a05", "h06_nu0_a02"]
#directories = ["h04_nu7_a100", "h04_nu7_a50", "h04_nu7_a167", "h04_nu7_a05"]
directories = ["h06_nu7_a05", "h06_nu7_a02"]

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = ".",
                         help = 'save directory (default: .)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = None,
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'maximum density (default: 1.1 times the max)')

    parser.add_argument('--not_log', dest = "arithmetic", action = 'store_true', default = False,
                         help = 'log or arithmetic (default: log)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 16,
                         help = 'fontsize of plot annotations (default: 16)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'fontsize of plot annotations (default: 3)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

def make_plot():
    # Data
    for directory in directories:
        xs = pickle.load(open("../%s/excess_mass_frames.p" % directory, "rb"))
        ys = pickle.load(open("../%s/excess_mass_values.p" % directory, "rb"))
        plot.plot(xs, ys, linewidth = args.linewidth, label = directory)

    # Annotate
    #title = r"$m_p = " + str(mass) + r" $ $M_J$, $\nu_{disk} = 10^{" + str(viscosity) + r"}$"
    plot.xlabel("Number of Planet Orbits", fontsize = args.fontsize)
    plot.ylabel("Excess Mass", fontsize = args.fontsize)
    #plot.title(title, fontsize = fontsize + 2)

    plot.legend(loc = "lower right")

    # Axes
    plot.xlim(args.r_lim[0], args.r_lim[1])
    #plot.ylim(0, 2)

    #plot.xscale("log")
    if args.arithmetic:
        pass
    else:
        plot.yscale("log")

    # Save + Close
    plot.savefig("excessMassesOverTime_set-%d.png" % directories[0])
    plot.show()

    plot.close()


make_plot()