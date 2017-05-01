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
start = 0
end = max_frame

full_range = np.arange(start, end)

save_rate = 50
sub_range = np.arange(start, end, save_rate)
#sub_range_two = np.arange(start - 1, end - 1, save_rate)

target_range = np.array([x for x in full_range if x not in sub_range])
#target_range = np.array([x for x in full_range if x not in sub_range and x not in sub_range_two])
pickle.dump(target_range, open("target_range_diagnostics.p", "wb"))

for prefix, new_prefix in zip(prefixes, new_prefixes):
    for frame in target_range:
        # Rename file
        old = "%s%d.dat" % (prefix, frame)
        new = "%s%d.dat" %(new_prefix, frame)

        if os.path.exists(old):
            shutil.move(old, new)
            #print old, new

