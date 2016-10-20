"""
retrieves active files for init, eos, and definitions for pluto

Usage:
python retrieveActiveFiles.py name
python retrieveActiveFiles.py name [ide]
"""

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

