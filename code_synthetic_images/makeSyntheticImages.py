"""
Generate synthetic images

Does
(1) Create tmp%06d directory
(2) Copy necessary files to this directory (use symbolic links)
(3) Name files appropriately
(4) Execute code
(5) Copy output back to original (or target) directory
(6) Delete tmp%06d directory
"""

import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import random

### Input Parameters ###

def new_argument_parser(description = "Generate input for synthetic images."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Synthetic Image Parameters
    parser.add_argument('-w', dest = "wavelength", type = float, default = 870,
                         help = 'wavelength (in um) (default: None)')

    parser.add_argument('-i', dest = "inclination", type = int, default = 0,
                         help = 'inclination (in degrees) (default: 0)')
    parser.add_argument('-T', dest = "temperature", type = int, default = 5770,
                         help = 'stellar temperature (in K) (default: 5770)')
    parser.add_argument('-R', dest = "radius", type = float, default = 6.96e10,
                         help = 'stellar radius (in cm) (default: 6.96e10)')
    parser.add_argument('-d', dest = "distance", type = int, default = 140,
                         help = 'distance (in pc) (default: 140)')

    # Save Parameters
    parser.add_argument('--id', dest = "id_number", type = int, default = 0,
                         help = 'id number (up to 4 digits) for this set of plot parameters (default: None)')
    parser.add_argument('--save', dest = "save_directory", default = None,
                         help = 'save directory (default: "lambda%d % wavelength")')

    return parser


###############################################################################

### Task Functions ###

def write_parameters():
	""" Step 0: Write synthetic image parameters to file """

	# Five Parameters
	stellar_radius = args.radius # in cm
	stellar_temperature = args.temperature # in K
	wavelength = args.wavelength # in um
	inclination = args.inclination
	distance = args.distance # in pc

	# Write to File
	txt = "%.1e\n%d\n%.1f\n%1.f\n%d" % (stellar_radius, stellar_temperature, wavelength, inclination, distance)
	fn = "parameters.dat"
	with open(fn, "w") as f:
		f.write(txt)

def setup_tmp_directory():
	""" Step 1: Make tmp directory """
	cwd = os.getcwd()

	# Make Directory
	random_number = random.randint(0, 999999)
	tmp_dir = "tmp%06d" % random_number
	os.mkdir(tmp_dir)

	# Fill it with necessary files
	necessary_files = []
	os.chdir(tmp_dir)
	os.symlink(target, name)

	os.chdir(cwd)


def full_procedure(frame):
    """ Every Step """
    density, sizes = retrieve_density(frame, directories)

    density = convert_units(density)
    density, sizes = polish(density, sizes)
    density = center_vortex(density)
    new_rad, new_theta, density = resample(density, new_num_rad = new_num_rad, new_num_theta = new_num_theta)
    density = interpolate_density(density, num_grains)
    output_density_txt(density, frame)
    output_density_pickle(density, frame)

    if frame == frame_range[0]:
        generate_secondary_files(new_rad, new_theta, sizes)
        save_id_parameters()

### Make Synthetic Images ###

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