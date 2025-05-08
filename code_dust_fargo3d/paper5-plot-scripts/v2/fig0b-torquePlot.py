"""
plot azimuthally averaged density
then, makes movies

Usage:
python plotAveragedDensity.py frame_number <== plot one frame
python plotAveragedDensity.py -m <== make a movie instead of plots (plots already exist)
python plotAveragedDensity.py -1 <<<===== Plots a sample
python plotAveragedDensity.py <== plot all frames and make a movie

"""

import sys, os, subprocess, csv
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

#import util
#import utilVorticity
#import azimuthal as az
#from readTitle import readTitle

#from advanced import Parameters
#from reader import Fields

#from colormaps import cmaps
#for key in cmaps:
#    plot.register_cmap(name = key, cmap = cmaps[key])

###############################################################################

### Input Parameters ###

def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--range', dest = "r_lim", type = float, nargs = 2, default = [0.3, 4.0],
                         help = 'radial range in plot (default: [r_min, r_max])')
    parser.add_argument('--max_y', dest = "max_y", type = float, default = 3.0,
                         help = 'maximum density (default: 1.1 times the max)')
    parser.add_argument('--y2_range', dest = "y2_range", type = float, nargs = 2, default = [0, 0.025],
                         help = 'range in y-axis (default: [-0.2, 0.2])')

    parser.add_argument('--log', dest = "log", action = 'store_true', default = False,
                         help = 'plot log y-axis (default: do not do it!)')

    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 21,
                         help = 'fontsize of plot annotations (default: 21)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 4,
                         help = 'fontsize of plot annotations (default: 3)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
#surface_density_zero = p.sigma0
#surface_density_power = p.sigmaslope

#taper_time = p.masstaper

scale_height = 0.1 #p.aspectratio
viscosity = 1.0e-5 #p.nu

#dt = p.ninterm * p.dt

#jupiter_mass = 1e-3
planet_mass = 3.0e-3 # fargo_par["PlanetMass"] / jupiter_mass
#accretion = fargo_par["Accretion"]
#taper_time = p.masstaper

### Get Input Parameters ###

save_directory = "."

# Plot Parameters (variable)
show = args.show

version = args.version
if args.r_lim is None:
    x_min = r_min; x_max = r_max
else:
    x_min = args.r_lim[0]; x_max = args.r_lim[1]
max_y = args.max_y
y2_range = args.y2_range

# Plot Parameters (constant)
fontsize = args.fontsize
linewidth = args.linewidth
dpi = args.dpi

###############################################################################

r_p = 1.0
r_hill = r_p * np.power(1.0 / 3.0 * planet_mass, 1.0 / 3.0)

# Helper Function
def my_model_J(r):

    distance = r - r_p

    #a = 0.41 * np.sqrt(6.0)
    a = 0.5
    torque = a * np.power(r, -0.5)

    if (np.abs(distance) <= 2.0 * r_hill):
      dr = np.abs(distance) / r_hill

      K = 5
      t = 1.5
      w = 0.15
      enhancement = 1.0 + 0.5 * (K - 1.0) * (1.0 - np.tanh((dr - (t + w) + 0.5 * w) / w))

      torque = enhancement * a

    return torque

###############################################################################

##### PLOTTING #####

alpha = 0.7

labelsize = 18
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    # Data (Our Model)
    rad = np.linspace(0.3, 4.0, 2000)
    torqueJs = np.array([my_model_J(r) for r in rad])

    x = rad
    y = torqueJs
    plot.plot(x, y, linewidth = linewidth, linestyle = "-", c = "b", zorder = 99, alpha = 0.8, label = "Model (This Work)")

    # Reference
    colors = ['k']
    num_jup = 3
    cases = ["torqueJ-Mt3"]
    for i, case_i in enumerate(cases):
       ref_density = []
       with open("%s.csv" % case_i, "r") as f:
          reader = csv.reader(f)
          for row in reader:
             ref_density.append(row)
       ref_density = np.array(ref_density).astype(np.float)

       #print np.shape(ref_density), ref_density

       xref = ref_density[:, 0] / 6.0 # Rp = 6 in Aoyama+Bai 23
       yref = ref_density[:, 1]
       #plot.plot(xref, yref, linewidth = linewidth, linestyle = "-", c = colors[i % 2], zorder = 2, alpha = 0.8, label = "A&B 2023")

    # Actual MHD
    colors = ['k', 'grey']
    #times = [20, 30, 50, 100, 140]
    times = []
    for i, time_i in enumerate(times):
       ref_density = []
       with open("Mt3Am3-t%d.csv" % time_i, "r") as f:
          reader = csv.reader(f)
          for row in reader:
             ref_density.append(row)
       ref_density = np.array(ref_density).astype(np.float)

       #print np.shape(ref_density), ref_density

       xref = ref_density[:, 0] / 6.0 # Rp = 6 in Aoyama+Bai 23
       yref = ref_density[:, 1]
       plot.plot(xref, yref, linewidth = linewidth, c = colors[i % 2], zorder = 2, label = "A&B 2023")

    # Reference
    ys = [-10, 10]

    xs1 = [r_p - r_hill, r_p - r_hill]
    xs2 = [r_p + r_hill, r_p + r_hill]
    plot.plot(xs1, ys, c = 'grey', linewidth = 2, linestyle = '--', alpha = 0.6)
    plot.plot(xs2, ys, c = 'grey', linewidth = 2, linestyle = '--', alpha = 0.6)

    xs3 = [r_p - 2.0 * r_hill, r_p - 2.0 * r_hill]
    xs4 = [r_p + 2.0 * r_hill, r_p + 2.0 * r_hill]
    plot.plot(xs3, ys, c = 'grey', linewidth = 1, linestyle = '--', alpha = 0.35)
    plot.plot(xs4, ys, c = 'grey', linewidth = 1, linestyle = '--', alpha = 0.35)

    #plot.legend(loc = "lower right", handlelength = 2.3, fontsize = fontsize - 6)

    # Axes
    min_y = 0.0
    if args.max_y is None:
        x_min_i = np.searchsorted(x, x_min)
        x_max_i = np.searchsorted(x, x_max)
        max_y = 1.1 * max(y[x_min_i : x_max_i])
    else:
        max_y = args.max_y

    #ax.set_xlim(x_min, x_max)
    ax.set_xlim(rad[0], rad[-1])
    #ax.set_xlim(0.5, 1.5)
    ax.set_ylim(min_y, max_y)

    planet_location = r_p
    plot.scatter([planet_location], [min_y + 1.05 * 10**(-3)], c = 'k', s = 75, alpha = 0.8, clip_on = False)

    # Annotate Axes
    #orbit = (dt / (2 * np.pi)) * frame

    unit = "r_\mathrm{p}"
    ax.set_xlabel(r"Radius [$%s$]" % unit, fontsize = fontsize)
    ax.set_ylabel(r"$dJ_\mathrm{MHD}/dr$ [$10^{-3}$ $\Sigma_0 r_\mathrm{p}^2 v_\mathrm{p} \Omega_\mathrm{p}$]", fontsize = fontsize)

    #if title is None:
    #    plot.title("Dust Density Map\n(t = %.1f)" % (orbit), fontsize = fontsize + 1)
    #else:
    #    plot.title("Dust Density Map\n%s\n(t = %.1f)" % (title, orbit), fontsize = fontsize + 1)

    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    y_text = 1.14

    if args.log:
        y_text = 2.15

    #title1 = r"$T_\mathrm{growth} = %d$ $\mathrm{orbits}$" % (taper_time)
    #title1 = r"$\Sigma_0 = %.3e$  $M_c = %.2f\ M_J$  $A = %.2f$" % (surface_density_zero, planet_mass, accretion)
    #title1 = r"$h/r = %.2f$     $\alpha \approx %s \times 10^{%d}$    $A = %.2f$" % (scale_height, alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2, accretion)
    #title2 = r"$t = %d$ $\mathrm{orbits}}$  [$M_\mathrm{p}\ =\ %.2f$ $M_\mathrm{Jup}$]" % (orbit, current_mass)
    title2 = "Torque Density"
    plot.title("%s" % (title2), y = 1.015, fontsize = fontsize + 1)
    #ax.text(x_mid, y_text * plot.ylim()[-1], title1, horizontalalignment = 'center', bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 2)

    # Text
    #text_mass = r"$M_\mathrm{p} = %d$ $M_\mathrm{Jup}$" % (int(planet_mass))
    #text_visc = r"$\alpha_\mathrm{disk} = 3 \times 10^{%d}$" % (int(np.log(viscosity) / np.log(10)) + 2)
    #text_comp = "Four-component"
    #text_comps = "(Torque, Mass Loss, Viscosity,"
    #text_comps2 = "& Extra Gap Torque)"
    #plot.text(0.95 * x_range + x_min, 0.92 * plot.ylim()[-1], text_comp, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.95 * x_range + x_min, 0.86 * plot.ylim()[-1], text_comps, fontsize = fontsize - 6, color = 'black', horizontalalignment = 'right')
    #plot.text(0.95 * x_range + x_min, 0.80 * plot.ylim()[-1], text_comps2, fontsize = fontsize - 6, color = 'black', horizontalalignment = 'right')
    #plot.text(-0.9 * box_size, 2, text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'left', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(0.9 * box_size, 2, text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))
    #plot.text(-0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_mass, fontsize = fontsize, color = 'black', horizontalalignment = 'right')
    #plot.text(0.84 * x_range / 2.0 + x_mid, y_text * plot.ylim()[-1], text_visc, fontsize = fontsize, color = 'black', horizontalalignment = 'left')


    # Inset
    x_in = [r_p - 3.0 * r_hill, r_p + 3.0 * r_hill]
    y_in = [min_y, max_y]
    #ax_inset = ax.inset_axes([0.5, 0.5, 0.47, 0.47], xlim = (x_in[0], x_in[1]), ylim = (y_in[0], y_in[1]))

    ax_inset = inset_axes(ax, width="45%", height="60%", borderpad = 2)

    plot.plot(x, y, linewidth = linewidth, linestyle = "-", c = "b", zorder = 99, alpha = 0.8, label = "Model (This Work)")

    plot.plot(xs1, ys, c = 'grey', linewidth = 2, linestyle = '--', alpha = 0.6)
    plot.plot(xs2, ys, c = 'grey', linewidth = 2, linestyle = '--', alpha = 0.6)

    plot.plot(xs3, ys, c = 'grey', linewidth = 1, linestyle = '--', alpha = 0.35)
    plot.plot(xs4, ys, c = 'grey', linewidth = 1, linestyle = '--', alpha = 0.35)

    plot.scatter([planet_location], [min_y + 1.05 * 10**(-3)], c = 'k', s = 75, alpha = 0.8, clip_on = False)

    ax_inset.set_xlim(x_in[0], x_in[1])
    ax_inset.set_ylim(min_y, max_y)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/torqueDensityProfile.png" % (save_directory)
    else:
        save_fn = "%s/v%04d_torqueDensityProfile.png" % (save_directory, version)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

##### Make Plot! #####

make_plot()
