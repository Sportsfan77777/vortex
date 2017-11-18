"""
label opacities with an id number

Usage: python labelOpacities.py id_number (previous_id)
"""

import os, shutil, glob, sys

def label_opacities(id_number, previous_id = None):
    """ add id_number to unlabeled (or designated) opacties"""

    # Find all opacity files
    default = "dustkappa_*.inp"
    if previous_id is not None:
        basename = "id%04d_%s" % (default, previous_id)
    else:
        basename = "%s" % default

    opacity_files = glob.glob(basename)

    # Re-label
    start = opacity_files[0].find("dustkappa")
    for opacity_file in opacity_files:
        src = opacity_file; dst = "id%04d_%s" % (id_number, opacity_file[start:])
        shutil.move(src, dst)


####### Label Opacities #######

if __name__ == "__main__":
    ### Main ###
    if len(sys.argv) == 2:
        # Usage: python labelOpacities.py id_number
        id_number = int(sys.argv[1]) # int
        label_opacities(id_number)
    elif len(sys.argv) == 3:
        # Usage: python labelOpacities.py id_number previous_id
        id_number = int(sys.argv[1]) # int
        previous_id = int(sys.argv[2]) # int
        label_opacities(id_number, previous_id = previous_id)
    else:
        print "Must supply 1 or 2 args\nUsage: python labelOpacities.py id_number (previous_id)"


