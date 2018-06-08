"""
analyze intensity
(1) print max intensity for range of files

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

from pylab import rcParams
from pylab import fromfile

import util
import square as sq
import azimuthal as az

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])

def new_argument_parser(description = "Plot convolved intensity maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('--ids', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')

    parser.add_argument('--ref', dest = "reference_id", type = int, default = None,
                         help = 'ref id number (up to 4 digits) to compare to (default: ids[0])')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Input Parameters ###
frames = util.get_frame_range(args.frames)
ids = util.get_frame_range(args.ids) # gets any range, really

reference_id = args.reference_id
if reference_id is None:
    reference_id = ids[0]

### Get ID%04d Parameters ###
fn = "id%04d_par.p" % reference_id
fargo_par = pickle.load(open(fn, "rb"))

num_rad = fargo_par["Nrad"]; num_theta = fargo_par["Nsec"]
r_min = fargo_par["Rmin"]; r_max = fargo_par["Rmax"]

jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
taper_time = fargo_par["MassTaper"]

surface_density_zero = fargo_par["Sigma0"] / 100
disk_mass = 2 * np.pi * surface_density_zero * (r_max - r_min) / jupiter_mass # M_{disk} = (2 \pi) * \Sigma_0 * r_p * (r_out - r_in)

scale_height = fargo_par["AspectRatio"]
viscosity = fargo_par["Viscosity"]

planet_radius = fargo_par["Radius"]

beam_size = fargo_par["Beam"]
wavelength = fargo_par["Wavelength"]
distance = fargo_par["Distance"]

arc_beam = beam_size * planet_radius / distance

###############################################################################

##### ANALYZE #####

for f, frame in enumerate(frames):
    print "t = %d" % frame

    reference_intensity_cart = util.read_data(frame, 'cartesian_intensity', fargo_par, id_number = reference_id)
    reference_max_intensity = np.max(reference_intensity_cart)

    for i, id_number in enumerate(ids):
        fn = "id%04d_par.p" % id_number
        fargo_par = pickle.load(open(fn, "rb"))

        intensity_cart = util.read_data(frame, 'cartesian_intensity', fargo_par, id_number = id_number)
        max_intensity = np.max(intensity_cart)
        percentage_loss = 100.0 * (1.0 - max_intensity / reference_max_intensity)

        print "id %04d: %.06f (%.1f%%)" % (id_number, max_intensity, percentage_loss)

    print






