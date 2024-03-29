"""
compare azimuthal intensity profiles at different beam sizes
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
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
    parser.add_argument('-b', dest = "beams", nargs = '+', type = int, default = [10, 20, 30, 40],
                         help = 'list of beams to compare (in AU) (default: [10, 20, 30, 40])')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalIntensityProfiles",
                         help = 'save directory (default: azimuthalIntensityProfiles)')

    # Old Format
    parser.add_argument('--old_res', dest = "old_res", type = int, nargs = 2, default = None,
                         help = 'resolution before re-sampling to calculate intensity (default: same)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", type = float, default = None,
                         help = 'max_y for each frame, or same for all (default: None)')
    parser.add_argument('--profiles', dest = "num_profiles", type = int, default = 1,
                         help = 'number of profiles (do not modify) (default: 1)')

    parser.add_argument('-n', dest = "normalize", action = 'store_false', default = True,
                         help = 'normalize by max (default: normalize)')

    parser.add_argument('-a', dest = "annotate", action = 'store_true', default = False,
                         help = 'annotate peak offsets at different thresholds (default: do not annotate)')
    
    # Plot Parameters (rarely need to change)
    parser.add_argument('--fontsize', dest = "fontsize", type = int, default = 19,
                         help = 'fontsize of plot annotations (default: 19)')
    parser.add_argument('--labelsize', dest = "labelsize", type = int, default = 16,
                         help = 'labelsize of plot annotations (default: 16)')
    parser.add_argument('--linewidth', dest = "linewidth", type = int, default = 3,
                         help = 'linewidths in plot (default: 3)')
    parser.add_argument('--alpha', dest = "alpha", type = float, default = 0.65,
                         help = 'line transparency in plot (default: 0.65)')
    parser.add_argument('--dpi', dest = "dpi", type = int, default = 100,
                         help = 'dpi of plot annotations (default: 100)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get ID%04d Parameters ###
fn = "id%04d_par.p" % (args.id_number)
fargo_par = pickle.load(open("../beam010/%s" % fn, "rb"))

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"]
dust_surface_density_zero = surface_density_zero / 100
disk_mass = 2 * np.pi * dust_surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

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

# Beam Choices
beams = args.beams

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Old Format
if args.old_res is None:
    old_num_rad = num_rad; old_num_theta = num_theta
else:
    old_num_rad = args.old_res[0]; old_num_theta = args.old_res[1]

# Plot Parameters (variable)
normalize = args.normalize

show = args.show
max_y = args.max_y
if normalize:
    max_y = 1

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

annotate = args.annotate

# Plot Parameters (constant)
fontsize = args.fontsize
labelsize = args.labelsize
linewidth = args.linewidth
alpha = args.alpha
dpi = args.dpi

rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

### Add new parameters to dictionary ###
fargo_par["rad"] = rad
fargo_par["theta"] = theta

###############################################################################

### Peak Offsets ###
frames = pickle.load(open("../beam010/id%04d_b10_t39_intensityFrames.p" % (id_number), 'rb'))
offsets1 = pickle.load(open("../beam010/id%04d_b10_t39_intensityPeaks.p" % (id_number), 'rb'))
#offsets2 = pickle.load(open("../beam020/id%04d_b20_t50_intensityPeaks.p" % (id_number), 'rb'))
#offsets3 = pickle.load(open("../beam030/id%04d_b30_t60_intensityPeaks.p" % (id_number), 'rb'))
#offsets4 = pickle.load(open("../beam040/id%04d_b40_t70_intensityPeaks.p" % (id_number), 'rb'))

###############################################################################

colors = ["#d11d1d", "#ef890b", "#4385ef", "#430aef", "k"]

alphas = [alpha, alpha, alpha, alpha, 1]
linestyles = ["--", "--", "-", "-", "-"]
linewidths = [linewidth, linewidth, linewidth + 1, linewidth + 1, linewidth + 1]

##### PLOTTING #####
def make_plot(frame, show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)

    ### Data ###
    for i, beam_i in enumerate(beams[::-1]):
        arc_beam_i = beam_i / distance
        label_i = r"$%.02f^{\prime\prime} (%2d \mathrm{\ AU})$" % (arc_beam_i, beam_i)

        intensity_polar = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number, directory = "../beam%03d" % beam_i)
        if normalize:
            intensity_polar /= np.max(intensity_polar)
        azimuthal_radius, azimuthal_profile = az.get_profiles(intensity_polar, fargo_par, args, shift = None)
        if normalize:
            azimuthal_profile /= np.max(azimuthal_profile)

        x = theta * (180.0 / np.pi)
        plot.plot(x, azimuthal_profile, linewidth = linewidths[i], c = colors[i], alpha = alphas[i], linestyle = linestyles[i], label = label_i)

    # Mark Planet (get shift first)
    if taper_time > 99.9:
        shift = az.get_lookup_shift(frame, directory = "../../../cm-size")
    else:
        dust_density = util.read_data(frame, 'dust', fargo_par, directory = "../../../cm-size")
        shift = az.get_azimuthal_peak(dust_density, fargo_par)

    if shift is None:
        planet_loc = theta[0]
    else:
        if shift < -len(theta):
            shift += len(theta)
        planet_loc = theta[shift] * (180.0 / np.pi)
    plot.scatter(planet_loc, 0, c = "k", s = 150, marker = "D", zorder = 100) # planet

    # Axes
    min_x = 0; max_x = 360
    plot.xlim(min_x, max_x)
    angles = np.linspace(min_x, max_x, 7)
    plot.xticks(angles)

    if max_y is None:
        plot.ylim(0, plot.ylim()[-1]) # No Input
    else:
        plot.ylim(0, max_y) # Input

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    current_mass = util.get_current_mass(orbit, taper_time, planet_mass = planet_mass)

    plot.xlabel(r"$\phi$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)
    plot.ylabel(r"$I$ / $I_\mathrm{max}$", fontsize = fontsize)

    # Title
    title = "\n" + r"$t$ $=$ $%.1f$   " % (orbit) + "[$m_p(t)$ $=$ $%.2f$ $M_J$]" % (current_mass)
    plot.title("%s" % (title), fontsize = fontsize + 1)

    # Legend
    #plot.legend(loc = "upper right", bbox_to_anchor = (1.38, 0.94)) # outside of plot
    if annotate:
        pass
    else:
        plot.legend(loc = "upper left")

    # Extra Annotation (Location, Legend Label)
    center_x = 1.38 * plot.xlim()[-1]
    top_y = plot.ylim()[-1]

    line = r"$\mathrm{Beam\ Diameters}$"
    #plot.text(center_x, 0.95 * top_y, line, fontsize = fontsize - 1, horizontalalignment = 'center')

    # Annotate Peak Offsets
    # if False:
    #     this_frame = np.searchsorted(frames, frame)
    #     offset1 = offsets1[this_frame]; offset2 = offsets2[this_frame]; offset3 = offsets3[this_frame]; offset4 = offsets4[this_frame]

    #     line4 = "b = 40: %.1f" % (offset4)
    #     line3 = "b = 30: %.1f" % (offset3)
    #     line2 = "b = 20: %.1f" % (offset2)
    #     line1 = "b = 10: %.1f" % (offset1)

    #     start_y = 0.08 * plot.ylim()[-1]; linebreak = 0.08 * plot.ylim()[-1]
    #     plot.text(180, start_y + 3.0 * linebreak, line4, fontsize = fontsize, horizontalalignment = 'center')
    #     plot.text(180, start_y + 2.0 * linebreak, line3, fontsize = fontsize, horizontalalignment = 'center')
    #     plot.text(180, start_y + 1.0 * linebreak, line2, fontsize = fontsize, horizontalalignment = 'center')
    #     plot.text(180, start_y + 0.0 * linebreak, line1, fontsize = fontsize, horizontalalignment = 'center')

    if annotate:
        this_frame = np.searchsorted(frames, frame - 0.1)
        offset0 = offsets1[this_frame - 3]; offset1 = offsets1[this_frame - 2]; offset2 = offsets1[this_frame - 1]; offset3 = offsets1[this_frame]; offset4 = offsets1[this_frame + 1]; offset5 = offsets1[this_frame + 2]

        line5 = "t = %d: %.1f (%.1f)" % (frame + 2, offset5, offset5 - offset4)
        line4 = "t = %d: %.1f (%.1f)" % (frame + 1, offset4, offset4 - offset3)
        line3 = "t = %d: %.1f (%.1f)" % (frame - 0, offset3, offset3 - offset2)
        line2 = "t = %d: %.1f (%.1f)" % (frame - 1, offset2, offset2 - offset1)
        line1 = "t = %d: %.1f (%.1f)" % (frame - 2, offset1, offset1 - offset0)

        start_y = 0.08 * plot.ylim()[-1]; linebreak = 0.08 * plot.ylim()[-1]
        plot.text(180, start_y + 4.0 * linebreak, line1, fontsize = fontsize, horizontalalignment = 'center')
        plot.text(180, start_y + 3.0 * linebreak, line2, fontsize = fontsize, horizontalalignment = 'center')
        plot.text(180, start_y + 2.0 * linebreak, line3, fontsize = fontsize, horizontalalignment = 'center')
        plot.text(180, start_y + 1.0 * linebreak, line4, fontsize = fontsize, horizontalalignment = 'center')
        plot.text(180, start_y + 0.0 * linebreak, line5, fontsize = fontsize, horizontalalignment = 'center')

    # Save, Show, and Close
    plot.tight_layout()

    png = "png"; pdf = "pdf"
    if version is None:
        save_fn = "%s/azimuthalIntensityProfiles_%04d.%s" % (save_directory, frame, png)
        pdf_save_fn = "%s/azimuthalIntensityProfiles_%04d.%s" % (save_directory, frame, pdf)
    else:
        save_fn = "%s/v%04d_azimuthalIntensityProfiles_%04d.%s" % (save_directory, version, frame, png)
        pdf_save_fn = "%s/v%04d_azimuthalIntensityProfiles_%04d.%s" % (save_directory, version, frame, pdf)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)
    plot.savefig(pdf_save_fn, bbox_inches = 'tight', format = "pdf")

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

# Iterate through frames

if len(frame_range) == 1:
    make_plot(frame_range[0], show = show)
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(make_plot, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            make_plot(frame)
