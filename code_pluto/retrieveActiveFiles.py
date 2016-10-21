"""
retrieves active files for init, eos, and definitions for pluto

Usage:
python retrieveActiveFiles.py name
python retrieveActiveFiles.py name [ide]
"""

import os
import shutil
import sys

# Pluto Code Directory
pluto_code = os.environ['PLUTO_DIR']

# Read Arguments
name = sys.argv[1]
if len(sys.argv) > 2:
    binary_options = sys.argv[2]
else:
    binary_options = "ide"

# Binary Options for Each File
move_init = "i" in binary_options
move_definitions = "d" in binary_options
move_eos = "e" in binary_options

# Do Not Use in Code Directory
current_directory = os.getcwd()
split = current_directory.split("/")
if split[-1] == "code_pluto":
    print "You are in the code directory. Do not use in this directory."
else:
    # Move Files
    if move_init:
        shutil.copy("%sinit_%s.c" % (pluto_code, name), "init_%s.c" % name)
    if move_definitions:
        shutil.copy("%sdefinitions_%s.h" % (pluto_code, name), "definitions_%s.h" % name)
    if move_eos:
        shutil.copy("%seos_%s.c" % (pluto_code, name), "eos_%s.c" % name)

