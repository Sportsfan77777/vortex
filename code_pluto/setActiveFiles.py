"""
updates current init, eos, and definitions files for pluto

Usage:
python setActiveFiles.py name
python setActiveFiles.py name [ide]
"""

import os
import shutil
import sys

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
        shutil.copy("init_%s.c" % name, "init.c")
    if move_definitions:
        shutil.copy("definitions_%s.h" % name, "definitions.h")
    if move_eos:
        shutil.copy("eos_%s.c" % name, "eos.c")