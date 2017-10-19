"""
Generate synthetic images

Does
(1) Create tmp%06d directory
[1a] Copy necessary files to this directory (use symbolic links)
[1b] Name files appropriately
(2) Execute code
[2a] Copy output to save directory
(2b) Delete tmp%06d directory
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

def setup_tmp_directory(frame):
	""" Step 1: Make tmp directory """
	cwd = os.getcwd()

	# Make Directory
	random_number = random.randint(0, 999999)
	tmp_dir = "tmp%04d" % frame
	os.mkdir(tmp_dir)

	# Files to move
	necessary_files = ["temperature.dat", "radial.dat", "grain.dat", "azimuthal.dat", "parameters.dat"]
	opacity_files = glob.glob("dustkappa_*.inp")
	density_file = "i%04d_gasddens%d.dat" % (id_number, frame)

	# Fill it with necessary files
	os.chdir(tmp_dir)
	for necessary_file in necessary_files:
		os.symlink("../%s" % necessary_file, necessary_file)

	# Fill it with opacity files	
	for opacity_file in opacity_files:
		os.symlink("../%s" % opacity_file, opacity_file)

	# Fill it with density file
	density_name = "sigmadust.dat"
	os.symlink("../%s" % density_file, density_name)

	os.chdir(cwd)
	return tmp_dir

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