"""
plot 2-D density maps

python plotDensityMaps.py
python plotDensityMaps.py frame_number
python plotDensityMaps.py -1 <<<===== Plots a sample
python plotDensityMaps.py -m
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
from multiprocessing import Array as mp_array
import argparse

import math
import numpy as np
from scipy.signal import find_peaks

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import gridspec

from pylab import rcParams
from pylab import fromfile

import util
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "extentsAndPeakCountsOverTime",
                         help = 'save directory (default: extentsOverTime)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", nargs = '+', type = float, default = None,
                         help = 'max_y for each frame, or same for all (default: None)')
    parser.add_argument('-s', dest = "sliver_width", type = float, default = 1.0,
                         help = 'number of scale heights in sliver (default: 1.0)')

    parser.add_argument('-t', dest = "threshold", type = float, default = 0.2,
                         help = 'threshold for measuring extent (default: 0.2)')

    parser.add_argument('--compare', dest = "compare", action = 'store_true', default = False,
                         help = 'compare the elongated vortex extents to the concentrated ones at the same threshold (default: do not compare)')

    parser.add_argument('--min', dest = "min_mass", type = float, default = 0.1,
                         help = 'minimum mass on plot (default: 0.1 Jupiter mass)')
    parser.add_argument('--max', dest = "max_mass", type = float, default = None,
                         help = 'maximum mass on plot (default: mass at last frame)')
    parser.add_argument('--delta', dest = "delta_mass", type = float, default = 0.1,
                         help = 'delta mass on plot (default: 0.1 Jupiter mass)')
    parser.add_argument('--minor', dest = "minor_delta_mass", type = float, default = None,
                         help = 'delta mass on plot (default: 0.1 Jupiter mass)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'labelsize of plot annotations (default: 16)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'linewidths in plot (default: 3)')
    parser.add_argument('--alpha', dest = "alpha", type = float, default = 0.35,
                         help = 'line transparency in plot (default: 0.35)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get ID%04d Parameters ###
fn = "id%04d_par.p" % args.id_number
fargo_par = pickle.load(open(fn, "rb"))

p = fargo_par["p"]

num_rad = p.ny; num_theta = p.nx
r_min = p.ymin; r_max = p.ymax

surface_density_zero = p.sigma0
dust_surface_density_zero = p.sigma0 * p.epsilon

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]

scale_height = p.aspectratio
viscosity = p.nu
taper_time = p.masstaper

dt = p.ninterm * p.dt

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

planet_radius = fargo_par["Radius"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
show = args.show
max_y = args.max_y

sliver_width = args.sliver_width

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

threshold = args.threshold

compare = args.compare

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

# Planet File
# Data
data = np.loadtxt("../../planet0.dat")
times = data[:, 0]; base_mass = data[:, 7] / jupiter_mass
accreted_mass = data[:, 8] / jupiter_mass

total_mass = base_mass + accreted_mass

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Data ###

def get_extent(args):
    # Extract args
    i, frame = args

    # Get data and measure extent
    intensity = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number)
    extent, azimuthal_profile = az.get_extent(intensity, fargo_par, normalize = True, threshold = threshold, sliver_width = sliver_width)

    # Count peaks
    peaks, _ = find_peaks(azimuthal_profile, height = threshold)
    peak_count = len(peaks)

    # Convert to degrees
    extent *= (180.0 / np.pi)

    print i, frame, extent, peak_count

    # Store
    extents[i] = extent
    peak_counts[i] = peak_count

###############################################################################

# Data
extents = mp_array("f", len(frame_range))
peak_counts = mp_array("f", len(frame_range))
pool_args = [(i, frame) for i, frame in enumerate(frame_range)]    

p = Pool(num_cores)
p.map(get_extent, pool_args)
p.terminate()


##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

colors = ['#d8db20', '#197229', '#519ba3', '#240f77'] # Ugly Yellow, Green, Slate Blue, Dark Blue

colors = ['#1f77b4', '#ff7f0e', '#be52e5', '#2ca02c'] # Blue, Orange, Purple, Green

size = 50
alpha_dots = 0.7

def make_plot(show = False):
    fig = plot.figure(figsize = (9, 6), dpi = dpi)
    #gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [7, 1.5], width_ratios = [7, 1], wspace = 0, figure = fig)
    gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [7, 1.5], width_ratios = [7, 1], figure = fig)
    ax = fig.add_subplot(gs[0, 0])

    # Plot
    x = frame_range
    y = np.array(extents)

    kernel = 5
    smooth_y = util.smooth(y, kernel)

    plot.plot(x, y, c = colors[1], linewidth = linewidth, alpha = alpha)
    plot.plot(x, smooth_y, c = colors[1], linewidth = linewidth)

    # Axes
    plot.xlim(x[0], x[-1])
    #ax.set_xticklabels([])

    angles = np.linspace(0, 360, 7)
    plot.yticks(angles)
    plot.ylim(0, 360)

    # Annotate Axes
    plot.ylabel(r"Azimuthal Extent $\Delta \phi$ $\mathrm{(degrees)}$", labelpad = 16, fontsize = fontsize + 2)

    threshold_text = r"$\frac{I_\mathrm{cut}}{I_\mathrm{max}}=%.2f$" % threshold
    #plot.text(0.98 * (x[-1] - x[0]) + x[0], 0.9 * plot.ylim()[-1], threshold_text, horizontalalignment = 'right', fontsize = fontsize - 4)
    plot.text(0.02 * (x[-1] - x[0]) + x[0], 0.05 * plot.ylim()[-1], threshold_text, horizontalalignment = 'left', fontsize = fontsize - 4)

    # Add text saying "More Elongated (up arrow)" and "More Compact (down arrow)"
    elongated_text = r"$\uparrow$More Elongated$\uparrow$"
    compact_text = r"$\downarrow$More Compact$\downarrow$"
    plot.text(0.65 * (x[-1] - x[0]) + x[0], 0.90 * plot.ylim()[-1], elongated_text, horizontalalignment = 'center', fontsize = fontsize - 3)
    plot.text(0.65 * (x[-1] - x[0]) + x[0], 0.05 * plot.ylim()[-1], compact_text, horizontalalignment = 'center', fontsize = fontsize - 3)

    #plot.legend(loc = "upper right", bbox_to_anchor = (1.28, 1.0)) # outside of plot
    #plot.legend(loc = "upper left") # outside of plot

    # Title
    x_min = plot.xlim()[0]; x_max = plot.xlim()[-1]
    x_range = x_max - x_min; x_mid = x_min + x_range / 2.0
    x_shift = 0.20; extra = 0.17
    y_text = 1.20; y_shift = 0.08

    alpha_coefficent = "3"
    if scale_height == 0.08:
        alpha_coefficent = "1.5"
    elif scale_height == 0.04:
        alpha_coefficent = "6"

    text1 = r"$h = %.2f$" % (scale_height)
    plot.text(x_min - x_shift * x_range, (y_text + y_shift) * plot.ylim()[-1], text1, horizontalalignment = 'left', fontsize = fontsize - 2)
    text2 = r"$\alpha \approx %s \times 10^{%d}$" % (alpha_coefficent, int(np.log(viscosity) / np.log(10)) + 2)
    plot.text(x_min - x_shift * x_range, (y_text) * plot.ylim()[-1], text2, horizontalalignment = 'left', fontsize = fontsize - 2)

    # Title #
    #title = r"$\mathrm{Azimuthal\ Extents}$"
    #title = r'$h = %.2f$   $\Sigma = %.3e$  (2-D)  [$%.3f^{\prime\prime}$]' % (scale_height, fargo_par["p"].sigma0, arc_beam)
    surface_density_base = 1.157e-4
    final_frame = 5000
    if final_frame > len(accreted_mass):
        final_frame = len(accreted_mass) - 1
    final_planet_mass = planet_mass + accreted_mass[final_frame]

    title1 = r'$%.3f^{\prime\prime}$' % (arc_beam)
    title2 = r"$\Sigma_0$ $/$ $\Sigma_\mathrm{base} = %.1f$    ($M_\mathrm{p} = %.2f$ $M_\mathrm{Jup}$)" % (surface_density_zero / surface_density_base, final_planet_mass)
    plot.title("[%s]\n%s" % (title1, title2), y = 1.30, fontsize = fontsize + 3, bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0))

    #### Peaks ####
    ax2 = fig.add_subplot(gs[1, 0])
    
    y2 = np.array(peak_counts)
    #plot.bar(x, y2, color = colors[2], edgecolor = colors[2], width = x[1] - x[0])
    plot.scatter(x, y2 - 0.15, color = colors[2], s = size, alpha = alpha_dots)

    # Axes
    plot.xlim(x[0], x[-1])

    counts = [0, 2, 4]
    plot.yticks(counts)
    plot.ylim(0, counts[-1])

    # Labels
    plot.xlabel(r"$t \mathrm{\ (planet\ orbits)}$", fontsize = fontsize + 2)
    plot.ylabel("# Peaks", fontsize = fontsize + 2)
    #plot.ylabel("# Peaks", fontsize = fontsize + 2, rotation = 270, labelpad = 25)

    #ax2.yaxis.set_label_position("right")
    #ax2.yaxis.tick_right()

    #### Histograms ####
    ax3 = fig.add_subplot(gs[0, 1])
    plot.hist(y, bins = np.linspace(0, 370, 38), cumulative = True, color = 'sienna', align = 'left', orientation = 'horizontal', histtype = 'stepfilled', density = True)

    ax3.set_xlim(0, 1)
    hist_ticks = np.linspace(0, 1, 11)
    hist_major_ticks = [0, 0.5, 1]
    ax3.set_xticks(hist_major_ticks)
    ax3.set_xticks(hist_ticks, minor = True)
    #ax3.set_xticklabels([])

    ax3.set_ylim(0, 360)
    ax3.set_yticks(angles)
    #ax3.set_yticklabels([])

    #ax3.set_title("Cumulative\nDistribution", fontsize = fontsize - 4)
    ax3.set_title("Probability\n" + r"($\leq \Delta \phi$)", fontsize = fontsize - 4)

    ax4 = fig.add_subplot(gs[1, 1])
    y2_adjusted = y2[:]
    y2_adjusted[y2 > 4] = 4
    plot.hist(y2_adjusted, bins = np.linspace(0, 4, 5) + 1e-8, cumulative = True, color = 'navy', orientation = 'horizontal', histtype = 'stepfilled', density = True)

    ax4.set_xlim(0, 1)
    ax4.set_xticks(hist_major_ticks)
    ax4.set_xticks(hist_ticks, minor = True)
    #ax4.set_xticklabels([])
    plot.xlabel("Cumulative\nFraction", fontsize = fontsize - 4)

    ax4.set_ylim(counts[0], counts[-1])
    ax4.set_yticks(counts)
    #ax4.set_yticklabels([])

    #### Add mass axis ####

    min_mass = args.min_mass; max_mass = args.max_mass; delta_mass = args.delta_mass
    if max_mass is None:
       max_mass = total_mass[frame_range[-1] - 1]

    mass_ticks = np.arange(min_mass, max_mass, delta_mass)

    def tick_function(masses):
        # For the secondary x-axis showing the planet mass over time
        tick_locations = np.zeros(len(masses))
        tick_labels = []

        for i, mass in enumerate(masses):
            #total_mass_jupiter = total_mass # in Jupiter masses
            times_i = az.my_searchsorted(total_mass, mass)

            #tick_times = times[times_i]

            print mass, times_i, len(times)

            tick_locations[i] = times[times_i]
            if delta_mass < 0.1:
                tick_labels.append("%.2f" % mass)
            else:
                tick_labels.append("%.1f" % mass)

        return tick_locations, tick_labels

    tick_locations, tick_labels = tick_function(mass_ticks)

    ax_twin = ax.twiny()
    ax_twin.set_xlim(ax.get_xlim())
    ax_twin.set_xticks(tick_locations)
    ax_twin.set_xticklabels(tick_labels)

    ax_twin.set_xlabel(r"$m_\mathrm{p}(t)$ [$M_\mathrm{J}$]", fontsize = fontsize - 1, labelpad = 8)

    if args.minor_delta_mass is not None:
        min_mass_minor = 0.1
        start_mass = total_mass[frame_range[0]]

        while min_mass_minor < start_mass:
            min_mass_minor += 0.05

        minor_mass_ticks = np.arange(min_mass_minor, max_mass, args.minor_delta_mass)
        minor_tick_locations, _ = tick_function(minor_mass_ticks)
        ax_twin.set_xticks(minor_tick_locations, minor = True)

    # Print counts
    print len(peak_counts)
    print (y2 == 1).sum(), (y2 == 2).sum(), (y2 == 3).sum(), (y2 >= 4).sum()

    num_one = 1.0 * (y2 == 1).sum() / len(peak_counts); num_two = 1.0 * (y2 == 2).sum() / len(peak_counts)
    num_three = 1.0 * (y2 == 3).sum() / len(peak_counts); num_four = 1.0 * (y2 >= 4).sum() / len(peak_counts)
    print "%.2f %.2f %.2f %.2f" % (num_one, num_two, num_three, num_four)

    # Save, Show, and Close
    current_directory = os.getcwd().split("/")[-3]
    current_beam = os.getcwd().split("/")[-1]
    if version is None:
        save_fn = "%s/extentsAndPeakCounts-%s-%s-%04d-%04d-%04d.png" % (save_directory, current_directory, current_beam, args.frames[0], args.frames[1], args.frames[2])
    else:
        save_fn = "%s/v%04d_extentsAndPeakCounts-%s-%s-%04d-%04d-%04d.png" % (save_directory, version, current_directory, current_beam, args.frames[0], args.frames[1], args.frames[2])
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi, pad_inches = 0.2)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)

    # Save Data
    ### TBD ###


##### Make Plots! #####

make_plot(show = show)

