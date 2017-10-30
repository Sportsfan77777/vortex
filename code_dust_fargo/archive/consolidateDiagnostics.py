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
prefixes = ['gasDiag1', 'gasDiag2', 'gasDiag3']
new_prefixes = ['rm_gasDiag1', 'rm_gasDiag2', 'rm_gasDiag3']

# Find number of frames
max_frame = util.find_max_frame()
num_frames = max_frame + 1

############# ############# ############# #############

# Set up range(s)
target_range = pickle.load(open("target_range_diagnostics.p", "rb"))

for new_prefix in new_prefixes:
    for frame in target_range:
        # Delete File
        old = "%s%d.dat" % (new_prefix, frame) # Note only re-named files can be deleted

        if os.path.exists(old):
            os.remove(old)
            #print old

# Delete file specifying target files
os.remove("target_range_diagnostics.p")