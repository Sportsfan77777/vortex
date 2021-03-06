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

prefixes = ['gasdens', 'gasvrad', 'gasvtheta']
new_prefixes = ['rm_gasdens', 'rm_gasvrad', 'rm_gasvtheta'] # re-named prefixes

density_files = glob.glob("*gasdens*.dat")

def find_max_frame():
    max_frame = 0
    for d_f in density_files:
        rename = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
        name = rename.split("_")[-1]
        frame_number = int(name[7:]) # just 999
        if frame_number > max_frame:
            max_frame = frame_number
    return max_frame

num_frames = find_max_frame() + 1

############# ############# ############# #############

rate = 5

for prefix in prefixes:
    for i in range(num_frames):
        if (i % rate) == 0:
            # All that is left should be re-numbered
            new_number = i / rate

            old = "%s%d.dat" % (prefix, i)
            new = "%s%d.dat" % (prefix, new_number)

            #shutil.move(old, new)
            print old, new


