import numpy as np
from matplotlib import pyplot as plot
from matplotlib import rcParams as rc
from matplotlib import ticker as ticker

import pickle

###############################################################################

master_directories = []

#master_directories.append(["h08_nu7_a167", "h08_nu7_a05", "h08_nu7_a02", "h08_nu7_a01"])
#master_directories.append(["h06_nu6_a50", "h06_nu6_a167", "h06_nu6_a05", "h06_nu6_a02"])
master_directories.append(["h06_nu7_a50", "h06_nu7_a167", "h06_nu7_a05", "h06_nu7_a02"])
#master_directories.append(["h06_nu0_a50", "h06_nu0_a167", "h06_nu0_a05", "h06_nu0_a02"])
#master_directories.append(["h04_nu7_a100", "h04_nu7_a50", "h04_nu7_a167", "h04_nu7_a05"])

###############################################################################

### Helper Functions ###

def get_final_mass(directory):
    planet_data = np.loadtxt("../%s/planet0.dat" % directory)
    times = planet_data[:, 6]
    masses = planet_data[:, 7]

    time_index = np.searchsorted(times, 3000)
    final_mass = masses[time_index]
    return final_mass

def record_lifetime(directory, cutoff = 0.2):
    # Helper
    def find_consecutive_ranges(array, values, cutoff, ranges = [], greater = True):
        if greater:
            test = lambda x : values[np.searchsorted(array, x)] > cutoff
        else:
            test = lambda x : values[np.searchsorted(array, x)] < cutoff

        for (match, group) in groupby(array, key = test):
            # Identify start and end of each range
            if match:
                start = next(group)
                end = start
                for end in group:
                    pass
                this_range = [start, end]
                ranges.append(this_range)
        return ranges

    frame_range = pickle.load(open("../%s/excess_mass_frames.p" % directory, "rb"))
    mass_over_time = pickle.load(open("../%s/excess_mass_values.p" % directory, "rb"))

    # Find Ranges When Vortex is Alive
    lifespans = find_consecutive_ranges(frame_range, mass_over_time, cutoff)

    # Add Up Lifetimes
    total_lifetime = 0
    for lifespan in lifespans:
        this_lifetime = lifespan[1] - lifespan[0]
        total_lifetime += this_lifetime

    # Print
    print "Lifespans: ", lifespans
    print "Total Lifetime: ", total_lifetime

    return lifetime


### PLOTTING ###
linewidth = 4
labelsize = 15
fontsize = 17
markersize = 9

colors = ["darkblue", "gold", "deepskyblue", "orange", "k"]
markers = ["s", "*", "p", "^", "h"]

def make_plot():
    # Data
    for i, directories in enumerate(master_directories):
        scale_height = float(directories[0].split("_")[0][1:])
        log_viscosity = float(directories[0].split("_")[1][1:])

        final_masses = []
        lifetimes = []

        for directory in directories:
            final_mass = get_final_mass(directory)
            lifetime = record_lifetime(directory, cutoff = 0.1)

            final_masses.append(final_mass)
            lifetimes.append(lifetime)

        xs = final_masses[::-1]
        ys = lifetimes[::-1]

        plot.plot(xs, ys, c = colors[i], marker = markers[i], markersize = markersize, linewidth = linewidth, label = case4_label)

    # Annotate
    plot.xlabel(r"Final Planet Mass (Jupiter masses)", fontsize = fontsize)
    plot.ylabel("Vortex Lifetime (planet orbits)", fontsize = fontsize)

    # Save + Show
    plot.savefig("present_lifetime_comparison.png", bbox_inches = "tight")
    plot.savefig("present_lifetime_comparison.pdf", bbox_inches = "tight", format = "pdf")
    plot.show()

make_plot()
