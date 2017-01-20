"""
save parameters from a single *.par file into a pickled dictionary
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
    elif line_sp[0][0] == "#":
        # Comment
        pass
    else:
        # Parse
        name = line_sp[0]
        entry = line_sp[1]

        # Store
        par[name] = entry

def store_planet_mass(line):
    """ stores only the planet's mass """
    line_sp = line.split()

    if len(line_sp) < 1:
        # Empty
        pass
    else:
        # Parse
        planet_name = line_sp[0]
        if planet_name == "Jupiter":
            par["PlanetMass"] = line_sp[2]

# Read file line by line
with open(par_file, "r") as f:
    for line in f:
        store(line)

### NEW: Add planet mass too ###
files = glob.glob("*.cfg") # usually Jup.cfg
cfg_file = files[0]

# Read file line by line
with open(cfg_file, "r") as f:
    for line in f:
        store_planet_mass(line) 

dict_name = "params.p"
p.dump(par, open(dict_name, "wb"))

