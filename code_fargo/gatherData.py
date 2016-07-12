"""
gets vortex start and end times for every case

Usage:
python gatherLifetimes.py
"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool
from multiprocessing import Array as mp_array

import math
import numpy as np
from scipy.ndimage import filters as ff

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot
from matplotlib import gridspec

from itertools import groupby

from pylab import rcParams # replace with rc ???
from pylab import fromfile

import util
from readTitle import readTitle

### Paths ###
base_path = "/rsgrps/kkratterstudents/mhammer/fargo_tests/fargo/mass_tapering_tests/%s/visc%d/evanescent/%s/"

def collect_data_for_case(case_number):
    if case_number == 1:
        # 1 MJ, visc6
        mass = "one_jupiter"
        visc = 6
        dir_name = "taper%d"
        tapers = [10, 250, 500]
    elif case_number == 2:
        # 1 MJ, visc7
        mass = "one_jupiter"
        visc = 7
        dir_name = "taper%d"
        tapers = [10, 500, 1000, 2000]
    elif case_number == 3:
        # 5 MJ, visc6
        mass = "five_jupiters"
        visc = 6
        dir_name = "taper%d_fold3and500"
        tapers = [10, 500, 1000, 2000]
    elif case_number == 4:
        # 5 MJ, visc7
        mass = "five_jupiters"
        visc = 7
        dir_name = "taper%d_fold3and500"
        tapers = [10, 1000, 2000, 4000]

    this_base_path = base_path % (mass, visc, dir_name)
    current_directory = os.getcwd() # store for later

    write_base = "case%s.p"

    lifetimes = []
    peak_densities = []
    trigger_masses = []
    for taper in tapers[1:]:
        # cd taper
        os.chdir(this_base_path % taper)

        # Read Lifetime
        #lifetime = pickle.load(open("lifetime.p", "rb"))
        #lifetimes.append(lifetime)

        # Read Peak Density
        #peak_density = pickle.load(open("peak_density.p", "rb"))
        #peak_densities.append(peak_density)

        # Read Trigger Mass
        factor = 1.0
        if case_number > 2:
            factor = 5.0

        trigger_dictionary = pickle.load(open("triggers.p", "rb"))
        trigger_time = trigger_dictionary[0.02]
        trigger_mass = factor * np.power(np.sin((np.pi / 2.0) * (1.0 * trigger_time / taper)), 2)

        trigger_masses.append(trigger_mass) # Use 0.02 as the trigger mass

        # cd back
        os.chdir(current_directory)

    # Write Tapers
    #taper_fn = write_base % "%d_tapers"
    #taper_fn = taper_fn % case_number
    #pickle.dump(tapers, open(taper_fn, "wb"))

    # Write Lifetime
    #lifetime_fn = write_base % "%d_lifetimes"
    #lifetime_fn = lifetime_fn % case_number
    #pickle.dump(lifetimes, open(lifetime_fn, "wb"))

    # Write Peak Density
    #peak_density_fn = write_base % "%d_densities"
    #peak_density_fn = peak_density_fn % case_number
    #pickle.dump(peak_densities, open(peak_density_fn, "wb"))

    # Write Trigger Mass
    trigger_fn = write_base % "%d_triggers"
    trigger_fn = trigger_fn % case_number
    pickle.dump(trigger_masses, open(trigger_fn, "wb"))

for i in range(4):
    collect_data_for_case(i + 1)