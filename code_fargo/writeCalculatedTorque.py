"""
write calculated torque to file
"""

import sys
import os
import subprocess
import pickle

import numpy as np

from pylab import fromfile

import util

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

max_frame = util.find_max_frame()
num_frames = max_frame + 1

#### Set Up and Binary Output (.npy) and Text Output (.dat) ####
## Binary ##
npy_fn = "calcTorque.npy"
npy_file = open(npy_fn, 'wb')

binary_array = np.zeros((8, num_frames)) - 1 # initialize to -1
np.save(npy_file, binary_array)
npy_file.close()

## Text ##
dat_fn = "calcTorque.dat"
dat_file = open(dat_fn, 'w')

# (1) Frame, (2) Total, (3) Inner, (4) Outer, 
# (5) InnerPositive, (6) InnerNegative, (7) OuterPositive, (8) OuterNegative
column_widths = 14 * np.ones(8, dtype = int)
column_widths[0] = 7

a = "Frame".center(column_widths[0])
b = "Total".center(column_widths[1])
c = "Inner".center(column_widths[2])
d = "Outer".center(column_widths[3])
e = "Inner (+)".center(column_widths[4])
f = "Inner (-)".center(column_widths[5])
g = "Outer (+)".center(column_widths[6])
h = "Outer (-)".center(column_widths[7])
first_line = "%s %s %s %s %s %s %s %s\n" % (a, b, c, d, e, f, g, h)
dat_file.write(first_line)
dat_file.close()

#### Calculate Data and Write to Files ####
for frame in range(num_frames):
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    torqueMap = util.torque(rad, theta, density)

    inner_torque_halves = util.inner_torque_contributions(rad, theta, torqueMap)
    outer_torque_halves = util.outer_torque_contributions(rad, theta, torqueMap)

    inner_torque = inner_torque_halves[0] + inner_torque_halves[1]
    outer_torque = outer_torque_halves[0] + outer_torque_halves[1]
    net_torque = inner_torque + outer_torque

    # Format into strings
    scaling = 10**6 # multiply by one million to make things readable
    a = ("%d" % frame).center(column_widths[0])
    b = ("%.8f" % (net_torque * scaling)).center(column_widths[1])
    c = ("%.8f" % (inner_torque * scaling)).center(column_widths[2])
    d = ("%.8f" % (outer_torque * scaling)).center(column_widths[3])
    e = ("%.8f" % (inner_torque_halves[0] * scaling)).center(column_widths[4])
    f = ("%.8f" % (inner_torque_halves[1] * scaling)).center(column_widths[5])
    g = ("%.8f" % (outer_torque_halves[0] * scaling)).center(column_widths[6])
    h = ("%.8f" % (outer_torque_halves[1] * scaling)).center(column_widths[7])

    line = "%s %s %s %s %s %s %s %s\n" % (a, b, c, d, e, f, g, h)

    # Write to File
    dat_file = open(dat_fn, 'a')
    dat_file.write(line)
    dat_file.close()

    # Fill in entries
    binary_array[0, frame] = frame
    binary_array[1, frame] = net_torque
    binary_array[2, frame] = inner_torque
    binary_array[3, frame] = outer_torque
    binary_array[4, frame] = inner_torque_halves[0]
    binary_array[5, frame] = inner_torque_halves[1]
    binary_array[6, frame] = outer_torque_halves[0]
    binary_array[7, frame] = outer_torque_halves[1]

    # Write to Binary
    npy_file = open(npy_fn, 'wb')
    np.save(npy_file, binary_array)
    npy_file.close()


