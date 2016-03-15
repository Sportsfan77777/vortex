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

import util

# All File Prefixes
prefixes = ['gasdens', 'gasvrad', 'gasvtheta']
new_prefixes = ['rm_gasdens', 'rm_gasvrad', 'rm_gasvtheta'] # re-named prefixes

fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    prefixes2D1D = [prefix + "1D" for prefix in prefixes]
    prefixes += prefixes2D1D

    new_prefixes2D1D = [new_prefix + "1D" for new_prefix in new_prefixes]
    new_prefixes += new_prefixes2D1D

# Find number of frames
max_frame = util.find_max_frame()
num_frames = max_frame + 1

############# ############# ############# #############

rate = 5

for new_prefix in new_prefixes:
    for i in range(num_frames):
        if (i % rate) == 0:
            # File should be kept
            pass
        else:
            # Delete File
            old = "%s%d.dat" % (new_prefix, i) # Note only re-named files can be deleted
        
            os.remove(old)
            #print old
