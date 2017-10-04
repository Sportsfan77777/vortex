"""
stores all function calls

Usage:
not to be called
"""

import os
import pickle

from pickleParameters import pickle_parameter_dictionary

#### Miscellaneous ####

def get_pickled_parameters():
    """ Retrieve parameters from *.par file """
    param_fn = "params.p"

    if not os.path.exists(param_fn):
        # create pickle if it doesn't exist
        pickle_parameter_dictionary()

    return pickle.load(open(param_fn, "rb"))

