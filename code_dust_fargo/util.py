"""
stores utility functions

Usage:
not to be called
"""

import os
import pickle

from pickleParameters import pickle_parameter_dictionary

#### Miscellaneous ####

def get_pickled_parameters(directory = "."):
    """ Retrieve parameters from *.par file """
    param_fn = "%s/params.p" % directory

    if not os.path.exists(param_fn):
        # create pickle if it doesn't exist
        pickle_parameter_dictionary(directory = directory)

    return pickle.load(open(param_fn, "rb"))

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
