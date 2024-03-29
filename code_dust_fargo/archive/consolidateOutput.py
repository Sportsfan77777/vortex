"""
reduces storage by keeping only 1/n of output files

(Step 1) python renameExtraOutput.py
(Step 2) python consolidateOutput.py
"""

import sys
import os
import shutil
import subprocess
import glob
import pickle
import numpy as np

import util

# All File Prefixes
prefixes = ['gasdens', 'gasvrad', 'gasvtheta', 'gasddens', 'gasdvrad', 'gasdvtheta', 'gasDiag1', 'gasDiag2', 'gasDiag3']
new_prefixes = ['rm_gasdens', 'rm_gasvrad', 'rm_gasvtheta', 'rm_gasddens', 'rm_gasdvrad', 'rm_gasdvtheta', 'rm_gasDiag1', 'rm_gasDiag2', 'rm_gasDiag3']

fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    prefixes2D1D = [prefix + "1D" for prefix in prefixes]
    prefixes2D1D_extra = ["gasdens.ascii_rad.", "gasvrad.ascii_rad.", "gasvtheta.ascii_rad.", \
                          "gaslabel", "gaslabel1D", "gaslabel.ascii_rad."]
    prefixes2D1D += prefixes2D1D_extra
    prefixes += prefixes2D1D

    new_prefixes2D1D = [new_prefix + "1D" for new_prefix in new_prefixes]
    new_prefixes2D1D_extra = ["rm_" + extra_prefix for extra_prefix in prefixes2D1D_extra]
    new_prefixes2D1D += new_prefixes2D1D_extra
    new_prefixes += new_prefixes2D1D

# Find number of frames
max_frame = util.find_max_frame()
num_frames = max_frame + 1

############# ############# ############# #############

# Set up range(s)
target_range = pickle.load(open("target_range.p", "rb"))

for new_prefix in new_prefixes:
    for frame in target_range:
        # Delete File
        old = "%s%d.dat" % (new_prefix, frame) # Note only re-named files can be deleted

        if os.path.exists(old):
            os.remove(old)
            #print old

# Delete file specifying target files
os.remove("target_range.p")