"""
stores utility functions

Usage:
not to be called
"""

import os
import pickle

import numpy as np

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

def read_data(frame, fn, directory = "."):
    """ read data"""
    # Dictionary
    basenames = {}
    basenames['gas'] = "gasdens%d.dat"; basenames['dust'] = "gasddens%d.dat"

    # Specific Data
    basename = basenames[fn] % frame
    density = (fromfile("%s/%s" % (directory, basename)).reshape(num_rad, num_theta)) * 100 # scale to gas density
    return density

def get_size_label(size):
    """ return label corresponding to size """
    sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001])
    size_labels = [r"$1$ $\rm{cm}$", r"$0.3$ $\rm{cm}$", r"$1$ $\rm{mm}$", r"$0.3$ $\rm{mm}$", r"$100$ $\rm{\mu m}$", r"$1$ $\rm{\mu m}$"]

    arg_size = np.abs(sizes - size).argmin() # find closest size
    return size_labels[arg_size]

def get_threshold(size):
    """ return label corresponding to size """
    ######## abstract this into a new method that takes an array as input!!!!!! ############
    sizes = np.array([1.0, 0.3, 0.1, 0.03, 0.01, 0.0001]) 
    thresholds = [5, 5, 2, 2, 2, 2]

    arg_size = np.abs(sizes - size).argmin() # find closest size
    return thresholds[arg_size]

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
