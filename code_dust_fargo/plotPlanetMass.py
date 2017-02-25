"""
plots planet mass over time
"""

import sys
import os
import subprocess
import pickle

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc

from scipy import signal as sig
from scipy.ndimage import filters as ff

from readTitle import readTitle

## Set file names ##
fargo_fn = "fargo2D1D"
if os.path.exists(fargo_fn):
    # fargo2D1D
    planet_fn = "bigplanet1.dat"
else:
    # fargo
    planet_fn = "bigplanet0.dat"

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

# Load Data (and choose subset) = x-axis
rate = 1 # If 1, choose all of the data. If >1, choose all_data / rate

data = np.loadtxt(planet_fn)
select = range(0, len(data[:,-1]), rate)
xs = (data[:, -4])[select] / (2 * np.pi) # Convert to num_orbits
masses = (data[:, 5])[select] # Planet Mass

# Plot Parameters
fontsize = 14
linewidth = 3

def make_plot():
    # Curves
    ys = masses / float(fargo_par["PlanetMass"])
    plot.plot(xs, ys, linewidth = linewidth)

    # Annotate
    this_title = readTitle()
    plot.title("%s" % this_title, fontsize = fontsize + 2)
    plot.xlabel("Orbit", fontsize = fontsize)
    plot.ylabel("Planet Mass / Initial Mass", fontsize = fontsize)

    # Limits
    #plot.ylim(1, 1.05)

    # Save and Close
    plot.savefig("planetMass.png", bbox_inches = 'tight')
    plot.show()

    plot.cla()


### PLOTTING ###

make_plot()