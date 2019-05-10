"""
stores utility functions

Usage:
not to be called
"""

import os
import pickle, glob

import numpy as np
from scipy.ndimage import filters as ff

from pylab import fromfile

from pickleParameters import pickle_parameter_dictionary

#### Miscellaneous ####

def get_pickled_parameters(directory = "."):
    """ Retrieve parameters from *.par file """
    param_fn = "%s/params.p" % directory

    if not os.path.exists(param_fn):
        # create pickle if it doesn't exist
        pickle_parameter_dictionary(directory = directory)

    return pickle.load(open(param_fn, "rb"))

def get_frame_range(frame_selection):
    """ return array of selected frames"""
    if len(frame_selection) == 1:
        frame_range = frame_selection
    elif len(frame_selection) == 2:
        start = frame_selection[0]; end = frame_selection[1]
        frame_range = range(start, end + 1)
    elif len(frame_selection) == 3:
        start = frame_selection[0]; end = frame_selection[1]; rate = frame_selection[2]
        frame_range = range(start, end + 1, rate)
    else:
        print "Error: Must supply 1, 2, or 3 frame arguments\nWith one argument, plots single frame\nWith two arguments, plots range(start, end + 1)\nWith three arguments, plots range(start, end + 1, rate)"
        exit()

    return frame_range

def read_data(frame, fn, fargo_par, id_number = None, version = None, directory = "."):
    """ read data"""
    ######## Get Parameters #########
    num_rad = fargo_par["Nrad"]
    num_theta = fargo_par["Nsec"]

    ########### Method ##############
    # Dictionary
    basenames = {}
    basenames['gas'] = "gasdens%d.dat"; basenames['dust'] = "gasddens%d.dat"; basenames['input_density'] = "id%04d_gasddens%d.p";
    basenames['intensity'] = "id%04d_intensity%04d.dat"; basenames['polar_intensity'] = "id%04d_intensityMap%04d.p"; basenames['cartesian_intensity'] = "id%04d_intensityCartGrid%04d.p"

    # Specific Data
    basename = basenames[fn]

    if "id" in basename:
        basename = basename % (id_number, frame)
    else:
        basename = basename % frame

    # Add Version
    if version is not None:
        basename = "v%04d_%s" % (version, basename)

    # Load properly based on extension
    ext = basename[basename.find("."):]
    if fn == 'intensity':
        data = (np.loadtxt("%s/%s" % (directory, basename))[:, -1]).reshape(num_rad, num_theta)
    elif ext == ".dat":
        data = (fromfile("%s/%s" % (directory, basename)).reshape(num_rad, num_theta))
    elif ext == ".p":
        data = pickle.load(open("%s/%s" % (directory, basename), "rb"))
    return data

def read_gas_data(frame, fargo_par, normalize = True, directory = "."):
    """ read dust data """
    ######## Get Parameters #########
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############
    data = read_data(frame, 'gas', fargo_par, directory = directory) 
    if normalize:
        data /= surface_density_zero
    return data

def read_dust_data(frame, fargo_par, normalize = True, directory = "."):
    """ read dust data """
    ######## Get Parameters #########
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############
    data = read_data(frame, 'dust', fargo_par, directory = directory) 
    if normalize:
        data /= (surface_density_zero / 100.0)
    return data

def read_merged_data(frame, num_cores, num_rad, num_theta, fn = 'gas', directory = "."):
    """ read data for outputs that were not merged """

    grid_name = "grid%03d.inf"

    # Dictionary
    basenames = {}
    basenames['gas'] = "gasdens%d_%d.dat"

    data = np.zeros((num_rad, num_theta))
    for i in range(num_cores):
        grid_name = grid_name % i
        basename = basenames[fn] % (frame, i)

        # Get size of each file
        with open(grid_name, 'r') as f:
            grid_lines = f.readlines()
            grid_numbers = grid_lines[1].split('\t')
            start_i = int(grid_numbers[1]); end_i = int(grid_numbers[2]) + 1; num_rad_i = end_i - start_i

        # Read data from each file
        data[start_i : end_i] = (fromfile(basename).reshape(num_rad_i, num_theta))

    return data

def save_merged_data(frame, num_cores, num_rad, num_theta, fn = 'gas', directory = "."):
    data = read_merged_data(frame, num_cores, num_rad, num_theta, fn = fn)

    # Dictionary
    basenames = {}
    basenames['gas'] = "gasdens%d.dat"

    basename = basenames[fn] % frame
    data.tofile(basename, format = "%.f")

    return data

def get_size(size_name):
    """ return number corresponding to size name """
    sizes = {}
    sizes["cm"] = 1.0; sizes["hcm"] = 0.3; sizes["mm"] = 0.1
    sizes["hmm"] = 0.03; sizes["hum"] = 0.01; sizes["um"] = 0.0001

    return sizes[size_name]

def get_size_name(size):
    """ return size name corresponding to size number """
    size_names = {}
    size_names[1.0] = "cm"; size_names[0.3] = "hcm"; size_names[0.1] = "mm"
    size_names[0.03] = "hmm"; size_names[0.01] = "hum"; size_names[0.0001] = "um"

    return size_names[size]

def get_size_label(size):
    """ return label corresponding to size """
    sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001]) ### <<<==== switch to dictionary!
    size_labels = [r"$1$ $\rm{cm}$", r"$3$ $\rm{mm}$", r"$1$ $\rm{mm}$", r"$0.3$ $\rm{mm}$", r"$100$ $\rm{\mu m}$", r"$1$ $\rm{\mu m}$"]

    arg_size = np.abs(sizes - size).argmin() # find closest size
    return size_labels[arg_size]

def get_stokes_number(size):
    """ return size name corresponding to size number """
    cm_ref = 0.076 # assumes r_p = 5 AU!
    return cm_ref * size

def get_threshold(size):
    """ return label corresponding to size """
    ######## abstract this into a new method that takes an array as input!!!!!! ############
    sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001]) 
    thresholds = [10.0, 5.0, 2.0, 1.7, 0.8, 0.7]

    arg_size = np.abs(sizes - size).argmin() # find closest size
    return thresholds[arg_size]

def get_current_mass(orbit, taper_time, planet_mass = 1.0):
    """ calculate mass at a particular orbit given a growth time up to a planet mass """
    if orbit >= taper_time:
        current_mass = planet_mass
    else:
        current_mass = np.power(np.sin((np.pi / 2) * (1.0 * orbit / taper_time)), 2) * planet_mass
    return current_mass

def find_max_frame():
    """ Get last orbit """
    fargo_fn = "fargo2D1D"
    if os.path.exists(fargo_fn):
        density_files = glob.glob("gasdens1D*.dat")
        number_index = 9 # the number starts at index 9 after 'gasdens1D'
    else:
        density_files = glob.glob("gasdens*.dat")
        number_index = 7 # the number starts at index 7 after 'gasdens'
    # Find the highest numbered file
    max_frame = 0
    for d_f in density_files:
        name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
        frame_number = int(name[number_index:]) # just 999
        if frame_number > max_frame:
            max_frame = frame_number
    return max_frame

def smooth(array, kernel_size):
    """ smoothing function """
    return ff.gaussian_filter(array, kernel_size)
