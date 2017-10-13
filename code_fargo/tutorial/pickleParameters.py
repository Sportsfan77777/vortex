"""
save parameters from a single *.par file into a pickled dictionary
"""

import pickle as p
import glob

def parse_parameters(par_dictionary):
    """ stores parameters in *.par file """

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

            # Convert to Proper Type (Integer, Float, or String?)
            is_int = False
            try:
                entry = int(entry) # Integer?
                is_int = True
            except ValueError:
                pass

            if not is_int:
                try:
                    entry = float(entry) # Float?
                except ValueError:
                    pass # String!

            # Store
            par_dictionary[name] = entry

    # Find *.par file (assume there is only one *.par file)
    files = glob.glob("*.par")
    par_file = files[0]

    # Read file line by line
    with open(par_file, "r") as f:
        for line in f:
            store(line)

def parse_planet_mass(par_dictionary):
    """ parse planet's mass """

    def store(line):
        """ stores only the planet's mass """
        line_sp = line.split()

        if len(line_sp) < 1:
            # Empty
            pass
        else:
            # Parse
            planet_name = line_sp[0]
            if planet_name == "Jupiter":
                par_dictionary["PlanetMass"] = float(line_sp[2])

    # Find *.cfg file (assume there is only one *.cfg file)
    files = glob.glob("*.cfg") # usually Jup.cfg
    cfg_file = files[0]

    # Read file line by line
    with open(cfg_file, "r") as f:
        for line in f:
            store(line)


def pickle_parameter_dictionary():
    """ dump all parameters into a pickle file """

    parameter_dictionary = {}

    parse_parameters(parameter_dictionary)
    parse_planet_mass(parameter_dictionary)

    dict_name = "params.p"
    p.dump(parameter_dictionary, open(dict_name, "wb"))

# Main
if __name__ == "__main__":
    pickle_parameter_dictionary()
