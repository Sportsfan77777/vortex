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
    parser.add_argument('frame', type = int,
                         help = 'select two frames to compare intensity (first is for T = 10, second is for T = 1000)')

    # Directory Selection
    parser.add_argument('--w1', dest = "wavelength1", type = float, default = 870,
                         help = 'wavelength (in um) (default: 870)')
    parser.add_argument('--w2', dest = "wavelength2", type = float, default = 3000,
                         help = 'wavelength (in um) (default: 300)')
    parser.add_argument('-b', dest = "beam_size", type = float, default = 10,
                         help = 'beam_size (in AU) (default: 10)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "azimuthalIntensityComparison",
                         help = 'save directory (default: azimuthalIntensityComparison)')

    # Plot Parameters (variable)
    parser.add_argument('--hide', dest = "show", action = 'store_false', default = True,
                         help = 'for single plot, do not display plot (default: display plot)')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of input parameters (default: 0)')
    parser.add_argument('-v', dest = "version", type = int, default = None,
                         help = 'version number (up to 4 digits) for this set of plot parameters (default: None)')

    parser.add_argument('--max_y', dest = "max_y", nargs = '+', type = float, default = None,
                         help = 'max_y for each frame, or same for all (default: None)')
    parser.add_argument('--profiles', dest = "num_profiles", type = int, default = 5,
                         help = 'number of profiles (default: 5)')
    parser.add_argument('-s', dest = "num_scale_heights", type = float, default = 8.0,
                         help = 'number of scale heights (default: 8.0)')

    parser.add_argument('-n', dest = "normalize", action = 'store_false', default = True,
                         help = 'normalize by max (default: normalize)')
    parser.add_argument('-t', dest = "threshold", type = float, default = None,
                         help = 'threshold for centering vortex with its center (default: varies with size)')
    
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
base_directory = "../taper1000/synthetic/lambda%04d/beam%03d"
default_directory = base_directory % (args.wavelength1, args.beam_size)

fn = "../%s/id%04d_par.p" % (default_directory, args.id_number)
fargo_par = pickle.load(open(fn, "rb"))

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
frame = args.frame

# Directories
directory1 = base_directory % (args.wavelength1, args.beam_size)
directory2 = base_directory % (args.wavelength2, args.beam_size)
directories = [directory1, directory2]

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

# Plot Parameters (variable)
normalize = args.normalize

show = args.show
max_y = args.max_y
if normalize:
    max_y = 1

num_profiles = args.num_profiles
num_scale_heights = args.num_scale_heights

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

id_number = args.id_number
version = args.version

threshold = args.threshold
#if threshold is None:
#    threshold = util.get_threshold(size)

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

### Helper Functions ###

def get_middle_profile(directory):
    # Change directories
    cwd = os.getcwd()
    os.chdir(directory)

    ### Data ###
    intensity_polar = util.read_data(frame, 'polar_intensity', fargo_par, id_number = id_number, directory = "lambda%04d/beam%03d" % (args.wavelength, args.beam_size))
    if normalize:
        intensity_polar /= np.max(intensity_polar)
    azimuthal_radii, azimuthal_profiles = az.get_profiles(intensity_polar, fargo_par, args, shift = None)

    # Middle Profile
    middle_i = (num_profiles - 1) / 2.0
    middle_profile = azimuthal_profiles[middle_i]

    # Return to previous directory
    os.chdir(cwd)

    return middle_profile

###############################################################################

##### PLOTTING #####

#colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
#          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
#          '#bcbd22', '#17becf']

#colors = ["#ff7f0e", "#9467bd", "#1f77b4", "#2ca02c", "#d62728"]
#dashes = [[3, 3], [42, 4], [10000, 1], [42, 4], [3, 3]]

#ns = [num_scale_heights / 2, num_scale_heights / 4, num_scale_heights / 4, num_scale_heights / 2]
#labels = [r"$\mathrm{-%.01f\ h}$" % ns[0], r"$\mathrm{-%.01f\ h}$" % ns[1], r"$r_\mathrm{c}$", r"$\mathrm{+%0.1f\ h}$" % ns[-2], r"$\mathrm{+%0.1f\ h}$" % ns[-1]]

colors = ["#1f77b4", "#9467bd"]
labels = [r"$%.02f$~\mathrm{mm}" % (wavelength1 / 100.0), r"$%.02f$~\mathrm{mm}" % (wavelength2 / 100.0)]

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)

    # Data
    middle_profiles = [get_middle_profile(directory_i) for directory_i in directories]

    # Plot
    x = theta * (180.0 / np.pi) - 180.0
    for i, middle_profile in enumerate(middle_profiles):
        plot.plot(x, middle_profile, linewidth = linewidth, c = colors[i], dashes = dashes[i], alpha = alpha, label = labels[i])

    # Save, Show, and Close
    png = "png"; pdf = "pdf"
    if version is None:
        save_fn = "%s/azimuthalIntensityBeamComparison_%04d.%s" % (save_directory, frame, png)
        pdf_save_fn = "%s/azimuthalIntensityBeamComparison_%04d.%s" % (save_directory, frame, pdf)
    else:
        save_fn = "%s/v%04d_azimuthalIntensityBeamComparison_%04d.%s" % (save_directory, version, frame, png)
        pdf_save_fn = "%s/v%04d_azimuthalIntensityBeamComparison_%04d.%s" % (save_directory, version, frame, pdf)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)
    plot.savefig(pdf_save_fn, bbox_inches = 'tight', format = "pdf")

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


def add_to_plot(frame, fig, ax, num_frames, frame_i):
    # Convert size to number
    size_name = "cm"
    size = util.get_size(size_name)

    # Taper
    if frame_i == 1:
        taper_time = 10
    else:
        taper_time = 1000

    ### Data ###
    middle_profiles = [get_middle_profile(directory_i) for directory_i in directories]

    # # Get Shift
    # dust_fargo_par = util.get_pickled_parameters(directory = "../cm-size") ## shorten name?
    # ######## Need to extract parameters, and add 'rad' and 'theta' ########
    # dust_rad = np.linspace(dust_fargo_par['Rmin'], dust_fargo_par['Rmax'], dust_fargo_par['Nrad'])
    # dust_theta = np.linspace(0, 2 * np.pi, dust_fargo_par['Nsec'])
    # dust_fargo_par['rad'] = dust_rad; dust_fargo_par['theta'] = dust_theta
    # gas_surface_density_zero = dust_fargo_par['Sigma0']

    # dust_density = util.read_data(frame, 'dust', dust_fargo_par, id_number = id_number, directory = "../cm-size")

    # # Shift gas density with center of dust density
    # shift = az.get_azimuthal_center(dust_density, dust_fargo_par, threshold = 10.0 * gas_surface_density_zero / 100.0)

    ### Plot ###
    # Profiles
    x = theta * (180.0 / np.pi) - 180.0
    for i, (radius, azimuthal_profile) in enumerate(zip(azimuthal_radii, azimuthal_profiles)):
        plot.plot(x, azimuthal_profile, linewidth = linewidth, c = colors[i], dashes = dashes[i], alpha = alpha, label = labels[i])

    # Mark Planet
    # if shift is None:
    #     planet_loc = theta[0]
    # else:
    #     if shift < -len(dust_theta):
    #         shift += len(dust_theta)
    #     planet_loc = dust_theta[shift] * (180.0 / np.pi) - 180.0
    #plot.scatter(planet_loc, 0, c = "k", s = 150, marker = "D", zorder = 100) # planet

    # Axes
    if taper_time < 10.1:
        # T = 10
        max_x = 180
    else:
        # T = 1000
        max_x = 180
    plot.xlim(-max_x, max_x)
    angles = np.linspace(-max_x, max_x, 7)
    plot.xticks(angles)

    if max_y is None:
        plot.ylim(0, plot.ylim()[-1]) # No Input
    else:
        plot.ylim(0, max_y) # Input

    # Annotate Axes
    time = fargo_par["Ninterm"] * fargo_par["DT"]
    orbit = (time / (2 * np.pi)) * frame
    current_mass = util.get_current_mass(orbit, taper_time, planet_mass = planet_mass)

    plot.xlabel(r"$\phi - \phi_\mathrm{center}$ $\mathrm{(degrees)}$", fontsize = fontsize + 2)

    if frame_i == 1:
        plot.ylabel(r"$I$ / $I_\mathrm{max}$", fontsize = fontsize)

    # Legend
    if frame_i == 2:
        plot.legend(loc = "upper right", bbox_to_anchor = (1.34, 0.94)) # outside of plot

    # Extra Annotation
    rc_line = r"$r_\mathrm{c} = %.02f$" % azimuthal_radii[(num_profiles - 1) / 2]
    plot.text(-170, 0.885 * plot.ylim()[-1], rc_line, fontsize = fontsize, horizontalalignment = 'left')

    if frame_i == 2:
        center_x = 1.38 * plot.xlim()[-1]
        top_y = plot.ylim()[-1]

        line = "Radii"
        #plot.text(center_x, 0.95 * top_y, line, fontsize = fontsize - 1, horizontalalignment = 'center')

    # Label
    #taper_title = r"$T_\mathrm{growth} = %d$" % taper_time
    #plot.text(0.95 * plot.xlim()[-1], 0.885 * plot.ylim()[-1], taper_title, fontsize = fontsize, color = 'black', horizontalalignment = 'right', bbox=dict(facecolor = 'white', edgecolor = 'black', pad = 10.0))

    # Title
    #title = "\n" + r"$t$ $=$ $%.1f$   " % (orbit) + "[$m_p(t)$ $=$ $%.2f$ $M_J$]" % (current_mass)
    #plot.title("%s" % (title), fontsize = fontsize + 1)

    # Return to previous directory
    os.chdir(cwd)

    
def old_make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)

    frame_str = ""
    for i, frame_i in enumerate(frame_range):
        add_to_plot(frame_i, fig, len(frame_range), i)
        frame_str += "%04d-" % frame_i
    frame_str = frame_str[:-1] # Trim last '_'

    #### Finish Plot ####
    #title = r"$%.03f^{\prime\prime} \times \ \ %.03f^{\prime\prime}$" % (arc_beam, arc_beam)
    #fig.suptitle(title, y = 0.97, verticalalignment = "bottom", bbox = dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad = 7.0), fontsize = fontsize + 4)

    # Save and Close
    plot.tight_layout()

    # Save, Show, and Close
    png = "png"; pdf = "pdf"
    if version is None:
        save_fn = "%s/azimuthalIntensityBeamComparison_%s.%s" % (save_directory, frame_str, png)
        pdf_save_fn = "%s/azimuthalIntensityTaperComparison_%s.%s" % (save_directory, frame_str, pdf)
    else:
        save_fn = "%s/v%04d_azimuthalIntensityBeamComparison_%s.%s" % (save_directory, version, frame_str, png)
        pdf_save_fn = "%s/v%04d_azimuthalIntensityBeamComparison_%s.%s" % (save_directory, version, frame_str, pdf)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)
    plot.savefig(pdf_save_fn, bbox_inches = 'tight', format = "pdf")

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


##### Make Plots! #####

make_plot(show = show)
