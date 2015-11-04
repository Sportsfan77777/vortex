"""
save parameters from *.par file into a pickled dictionary
"""

import pickle as p
import glob

# Find *.par file (assume there is only one *.par file)

par = {}

files = glob.glob("*.par")
par_file = files[0]

def store(line):
	""" For a line 'a b', stores 'b' under 'a' in the dictionary """
	line_sp = line.split()

	if len(line_sp) < 2:
		# Empty
		pass
	if line_sp[0][0] == "#":
		# Comment
		pass

	# Parse
	name = line_sp[0]
	entry = line_sp[1]

	# Store
	par[name] = entry

# Read file line by line