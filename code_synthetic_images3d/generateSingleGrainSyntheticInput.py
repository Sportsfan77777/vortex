"""
Usage:
python generateSyntheticInput.py

Does
(1) Options to carve out inner cavity (default: yes) and scale dust sizes (default: no)
(2) Centers vortices (using cm, hcm, and mm peaks as a guide -- For T = 1000, center assumed to be fixed)
(3) Re-samples at a lower resolution
(4) Interpolates dust to continuous distribution of grain sizes (e.g. 100)
(5) Convert density distributions from code units to cgs
(6) Generates secondary input: radial.dat, azimuthal.dat, temperature.dat, grain.dat
(7) Outputs "gasddens%d.dat", a text file input for the synthetic image generator (access with "ln -s gasddens%d.dat" density.dat)
(8) Output "gassdens%d.p", a pickle file that can be used to look at the overall dust density distribution (such as to look at the gas-to-dust ratio over time)
(9) Save extra parameters in new dictionary: "id%04d_par.p"

Assumes directory structure with adjacent access to directories containing different grain sizes
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np
from scipy import interpolate as sp_int

import argparse
from pylab import fromfile

import util
import azimuthal as az
from labelOpacities import label_opacities

from advanced import Parameters
from reader import Fields

size_names = ["hcm"] # , "hmm", hum", "um"]
sizes = np.array([0.3]) #, 0.03, 0.01, 0.0001])

### Input Parameters ###

def new_argument_parser(description = "Generate input for synthetic images."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Select Dust
    parser.add_argument('-n', dest = "n", type = float, default = 1,
                         help = 'choose grain number (default: 1)')
    parser.add_argument('--size', dest = "size", type = float, default = 0.3,
                         help = 'choose size (default: 0.3)')
    parser.add_argument('--name', dest = "name", default = "hcm",
                         help = 'choose size name (default: hcm)')

    # System Parameters
    parser.add_argument('-m', dest = "mass", type = float, default = 1.0,
                         help = 'mass of star in solar masses (default: 1.0)')
    parser.add_argument('-r', dest = "radius", type = float, default = 20.0,
                         help = 'radius of planet in AU (default: 20.0)')

    # Centering
    parser.add_argument('--shift', dest = "center", default = "away",
                         help = 'center options: threshold, cm-threshold, away (from minimum), cm-away, peak, cm-peak, lookup (default: lookup)')

    # Contribution Test
    parser.add_argument('-z', dest = "zero", type = int, default = -1,
                         help = 'index of grain to leave out for contribution test (default: -1)')

    # Interpolation
    parser.add_argument('--interpolate', dest = "interpolate", action = 'store_true', default = False,
                         help = 'interpolate or not (default: do not interpolate)')

    # Save Parameters
    parser.add_argument('-p', dest = "number_density_power", type = float, default = 3.5,
                         help = 'negative power in grain size power law (default: 3.5)')
    parser.add_argument('--scale-density', dest = "scale_density", type = float, default = 100,
                         help = 'scaling of grain size distribution (default: 1)')
    parser.add_argument('--scale-sizes', dest = "scale_sizes", type = float, default = 1,
                         help = 'scaling of grain size distribution (default: 1)')

    parser.add_argument('-s', dest = "new_res", nargs = 2, type = int, default = None,
                         help = 're-sample resolution (default: [Nrad, Nsec])')
    parser.add_argument('-t', dest = "new_range", nargs = 2, type = float, default = None,
                         help = 're-sample range (default: [Rmin, Rmax])')
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of plot parameters (default: None)')
    parser.add_argument('--save', dest = "save_directory", default = ".",
                         help = 'save directory (default: ".")')
    parser.add_argument('--separate', dest = "save_separate", action = 'store_true', default = False,
                         help = 'save density for num_grains grain sizes (default: False)')

    parser.add_argument('-o', dest = "make_opacities", action = 'store_true', default = False,
                         help = 'make opacities (default: False)')

    return parser

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get FARGO Parameters ###
p = Parameters()
fargo_par = util.get_pickled_parameters()

num_rad = p.ny; num_theta = p.nx
r_min = p.ymin; r_max = p.ymax

surface_density_zero = p.sigma0
dust_surface_density_zero = p.sigma0 * p.epsilon

planet_mass = 1.0
taper_time = p.masstaper

scale_height = p.aspectratio
viscosity = p.nu

dt = p.ninterm * p.dt

rad = np.linspace(r_min, r_max, num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)

### Get Input Parameters ###

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Dust
sizes[0] = args.size
size_names[0] = args.name

# System Parameters
mass = args.mass # (in solar masses)
radius = args.radius # radius of planet (in AU)

# Centering
center = args.center

# Interpolation
interpolate = args.interpolate

# Contribution Test
zero = args.zero

# Save Parameters
num_grains = 1
number_density_power = -args.number_density_power # Note: negative of input
scale_density = args.scale_density
scale_sizes = args.scale_sizes

if args.new_res is None:
    new_num_rad = num_rad; new_num_theta = num_theta
else:
    new_num_rad = args.new_res[0]; new_num_theta = args.new_res[1]
if args.new_range is None:
    new_r_min = r_min; new_r_max = r_max
else:
    new_r_min = args.new_range[0]; new_r_max = args.new_range[1]

id_number = args.id_number
save_directory = args.save_directory
save_separate = args.save_separate

make_opacities = args.make_opacities

### Miscellaneous ###
# Constants (G, \mu, m_p, k_b)
G = 6.67 * 10**-8; mu = 2.34; mp = 1.67 * 10**-24; kb = 1.38 * 10**-16

# Units
mass_unit = mass * (1.988425 * 10**33) # (solar mass / g)
radius_unit = radius * (1.496 * 10**13) # (AU / cm)
density_unit = mass_unit / radius_unit**2 # unit conversion factor

###############################################################################

### Add new parameters to dictionary ###
fargo_par["p"] = p

fargo_par["rad"] = rad
fargo_par["theta"] = theta

fargo_par["Nrad"] = num_rad; fargo_par["Nsec"] = num_theta

fargo_par["new_num_rad"] = new_num_rad; fargo_par["new_num_theta"] = new_num_theta
fargo_par["new_r_min"] = new_r_min; fargo_par["new_r_max"] = new_r_max

###############################################################################

### Task Functions ###

def retrieve_density(frame, size_names):
    """ Step 0: Retrieve density """
    density = np.zeros((num_rad, num_theta, len(sizes)))
    starting_sizes = sizes

    for i, size_name in enumerate(size_names):
        directory = "."

        if zero == i:
            # For debugging purposes to test contribution of a single grain size
            pass
        else:
            density[:, :, i] = util.read_dust_data(frame, fargo_par, normalize = False, directory = directory, n = args.n)

    return density, starting_sizes

def polish(density, sizes, cavity_cutoff = 0.92, scale_density = 1, scale_sizes = 1):
    """ Step 1: get rid of inner cavity and scale dust densities to different grain size """
    # Cavity
    if cavity_cutoff is not None:
        cavity_cutoff_i = np.searchsorted(rad, cavity_cutoff)
        density[:cavity_cutoff_i, ] /= 100.0

    # Scale
    density *= scale_density
    sizes *= scale_sizes

    return density, sizes

def center_vortex(density, frame, reference_density = None):
    """ Step 2: center the vortex so that the peak is at 180 degrees """
    if taper_time < 10.1:
        for i, size in enumerate(sizes):
            shift_i = az.get_azimuthal_peak(density[:, :, i], fargo_par)
            density[:, :, i] = np.roll(density[:, :, i], shift_i, axis = 1)
        return density

    elif taper_time > 99.9:
        for i, size_name in enumerate(size_names):
            if center == "lookup":
                pass
                #shift_i = az.get_lookup_shift(frame, directory = "../%s-size" % size_name)
            elif center == "threshold":
                threshold = util.get_threshold(size) * surface_density_zero
                shift_i = az.get_azimuthal_center(density[:, :, i], fargo_par, threshold = threshold)
            elif center == "away":
                shift_i = az.shift_away_from_minimum(reference_density, fargo_par)
            else:
                print "invalid centering option"
            density[:, :, i] = np.roll(density[:, :, i], shift_i, axis = 1)
        return density

def resample(density, new_num_rad = 400, new_num_theta = 400):
    """ Step 3: lower resolution (makes txt output smaller) """

    new_density = np.zeros((new_num_theta, new_num_rad, len(sizes)))

    new_rad = np.linspace(new_r_min, new_r_max, new_num_rad)
    new_theta = np.linspace(0, 2 * np.pi, new_num_theta)

    for i, _ in enumerate(sizes):
        interpolator = sp_int.interp2d(theta, rad, density[:, :, i]) # Careful: z is flattened!
        new_density[:, :, i] = (interpolator(new_theta, new_rad)).T # Note: result needs to be transposed
    
    return new_rad, new_theta, new_density

def interpolate_density(density, num_grains):
    """ Step 4 - Option A: interpolate to more grain sizes """

    ### New Grain Sizes and Ranges ###
    log_interpolated_sizes = np.linspace(np.log10(min(sizes)), np.log10(max(sizes)), num_grains)
    interpolated_sizes = np.power(10.0, log_interpolated_sizes) # to be used in size interpolation

    # Grain Ranges (for power-law distribution)
    da = log_interpolated_sizes[1] - log_interpolated_sizes[0]
    log_interpolated_ranges = np.linspace(np.log10(min(sizes)) - 0.5 * da, np.log10(max(sizes)) + 0.5 * da, num_grains + 1)
    interpolated_ranges = np.power(10.0, log_interpolated_ranges)

    ### Interpolate ### 
    size_interpolator = sp_int.interp1d(sizes, density, axis = -1)
    unweighted_interpolated_density = size_interpolator(interpolated_sizes)

    # Scale to a Power Law Distribution
    surface_density_power = number_density_power + 3.0

    def integrate_size_distribution(x, y):
        """ Determine weight of each dust size bin """
        def f(z):
            return z ** (surface_density_power + 1.0) # +1 for integration
        return f(y) - f(x)

    range_a = interpolated_ranges[:-1]
    range_b = interpolated_ranges[1:]
    size_weights = [integrate_size_distribution(x, y) for (x, y) in zip(range_a, range_b)]
    power_law = size_weights / np.sum(size_weights)

    interpolated_density = power_law[None, :] * unweighted_interpolated_density

    # Normalize to initial total gas mass

    return interpolated_density, interpolated_sizes

def distribute_density(density):
    """ Step 4 - Option B: distribute density according to a power law """
    # Reverse Sizes (may add different options later)
    new_sizes = sizes[::-1]

    # Reverse Density
    reversed_density = density[:, :, ::-1]

    # Grain Ranges (for power-law distribution)
    log_sizes = np.log10(sizes[::-1])
    range_a = log_sizes[:-1]
    range_b = log_sizes[1:]
    log_middle_midpoint_sizes = [np.mean([a, b]) for (a, b) in zip(range_a, range_b)]

    log_midpoint_sizes = np.zeros(len(log_middle_midpoint_sizes) + 2)
    log_midpoint_sizes[1:-1] = log_middle_midpoint_sizes
    log_midpoint_sizes[0] = log_sizes[0] # - 0.5 * (middle_midpoint_sizes[0] - log__sizes[0])
    log_midpoint_sizes[-1] = log_sizes[-1] #+ 0.5 * (middle_midpoint_sizes[-1] - log_sizes[-1])

    midpoint_sizes = np.power(10.0, log_midpoint_sizes)
    midpoint_sizes[0] = 0 # Fix smallest size to 0
    dust_bins = np.diff(midpoint_sizes)

    # Scale to a Power Law Distribution
    surface_density_power = number_density_power + 3.0

    def weighting_distribution(x):
        """ Determine weight of each dust size bin """
        return x ** (surface_density_power)

    size_weights = [weighting_distribution(s) * dust_bins[i] for (i, s) in enumerate(sizes[::-1])]
    power_law = size_weights / np.sum(size_weights)

    distributed_density = power_law[None, :] * reversed_density

    return distributed_density, new_sizes

def convert_units(density):
    """ Step 5: convert density from code units to cgs units """
    density *= density_unit
    return density

def generate_secondary_files(rad, theta, new_sizes):
    """ Step 6: write other *.dat files + make opacities """

    def temperature(R):
        """ r """
        R *= radius_unit # convert to cgs

        # (mu * mp) (h**2 \omega **2) / (k_b)
        omega_sq = (G * mass_unit) / (R**3)
        scale_height_cgs = scale_height * (R)
        return (mu * mp) * (scale_height_cgs**2 * omega_sq) / (kb) 
        
    temperatures = np.array([temperature(r) for r in rad])

    # Save Files
    np.savetxt("%s/id%04d_radial.dat" % (save_directory, id_number), rad * radius_unit)
    np.savetxt("%s/id%04d_azimuthal.dat" % (save_directory, id_number), theta)
    np.savetxt("%s/id%04d_grain.dat" % (save_directory, id_number), new_sizes)
    np.savetxt("%s/id%04d_temperature.dat" % (save_directory, id_number), temperatures)

    # Make and Label Opacities
    if make_opacities:
        np.savetxt("%s/grain.dat" % (save_directory), new_sizes) # need grain size distribution to make opacities
        np.savetxt("%s/temperature.dat" % (save_directory), temperatures) # not sure you need temperature

        os.system("./makeopac")
        label_opacities(id_number)

def output_density_txt(density, frame):
    """ Step 7: output txt file """
    interleaved_density = density.flatten('F') # interleave to 1-d

    fn = "%s/id%04d_gasddens%d.dat" % (save_directory, id_number, frame)
    np.savetxt(fn, interleaved_density)

def output_density_pickle(density, frame):
    """ Step 8: output pickle file """
    # Save Separated Density
    if save_separate:
        fn = "%s/id%04d_separated_gasddens%d.p" % (save_directory, id_number, frame)
        pickle.dump(density, open(fn, 'wb'))

    # Save Composite Density
    composite_density = np.sum(density, axis = -1)

    fn = "%s/id%04d_gasddens%d.p" % (save_directory, id_number, frame)
    pickle.dump(composite_density, open(fn, 'wb'))

def save_id_parameters():
    """ Step 9: Save parameters associated with this id number """
    id_par = fargo_par.copy()

    id_par["id"] = id_number

    id_par["Mass"] = mass # in solar masses
    id_par["Radius"] = radius # in AU
    id_par["MassUnit"] = mass_unit
    id_par["RadiusUnit"] = radius_unit

    id_par["Rmin"] = new_r_min
    id_par["Rmax"] = new_r_max
    id_par["Nrad"] = new_num_rad
    id_par["Nsec"] = new_num_theta
    id_par["Sigma0"] *= density_unit

    id_par["old_Nrad"] = num_rad
    id_par["old_Nsec"] = num_theta

    id_par["NumberDensityPower"] = number_density_power

    dict_name = "id%04d_par.p" % id_number
    pickle.dump(id_par, open(dict_name, "wb"))

def full_procedure(frame):
    """ Every Step """
    gas_density = util.read_gas_data(frame, fargo_par, normalize = False)
    density, sizes = retrieve_density(frame, size_names)

    density, sizes = polish(density, sizes, scale_density = scale_density, scale_sizes = scale_sizes)
    if center != "off":
        density = center_vortex(density, frame, reference_density = gas_density)
    new_rad, new_theta, density = resample(density, new_num_rad = new_num_rad, new_num_theta = new_num_theta)

    if interpolate:
        density, new_sizes = interpolate_density(density, num_grains)
    else:
        density, new_sizes = distribute_density(density)

    density = convert_units(density)
    output_density_txt(density, frame)
    output_density_pickle(density, frame)

    if frame == frame_range[0]:
        generate_secondary_files(new_rad, new_theta, new_sizes)
        save_id_parameters()

###############################################################################

### Generate Synthetic Input ###

# Iterate through frames

if len(frame_range) == 1:
    full_procedure(frame_range[0])
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(full_procedure, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            full_procedure(frame)


