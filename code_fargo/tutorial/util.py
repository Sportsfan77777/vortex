"""
stores all function calls

Usage:
not to be called
"""

import os
import subprocess
import glob
import time
import pickle

import numpy as np

#### Miscellaneous ####

def get_pickled_parameters():
    """ Retrieve parameters from *.par file """
    param_fn = "params.p"

    if not os.path.exists(param_fn):
        # create pickle if it doesn't exist
        command = "python pickleParameters.py"
        subprocess.Popen(command.split())
        time.sleep(1) # one sec -- give it time to execute

    return pickle.load(open(param_fn, "rb"))

