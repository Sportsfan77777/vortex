"""
save parameters from a single *.par ini into a pickled dictionary
"""

import pickle as p
import glob

# Find *.ini file (assume there is only one *.ini file)

par = {}

files = glob.glob("*.ini")
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
        entry = line_sp[1:]

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

dict_name = "params.p"
p.dump(par, open(dict_name, "wb"))

