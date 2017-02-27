"""
save pickle file of dictionary containing
(1) cube: cube of azimuthally averaged density over time
(2) scale_height: scale height
(3) rad: used_radii
(4) id: identifier dictionary --- just fargo_par

Usage:
python gatherAveragedDensityCube.py
** Then, compare with another cube. **
** Then, make a movie. **

"""

import sys
import os
import subprocess
import glob
import pickle
from multiprocessing import Pool

import math
import random
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import rcParams # replace with rc ???
from pylab import fromfile

import util
from readTitle import readTitle

### Get FARGO Parameters ###
# Create param file if it doesn't already exist
pickled = util.pickle_parameters()
param_fn = "params.p"
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])

#### Gather Data ####

## Use These Frames ##
rate = 5
start = 0
max_frame = 2000 #util.find_max_frame()
frame_range = np.array(range(start, max_frame + 1, rate))
num_frames = len(frame_range)

## Set Up Cube ##
averagedDensityCube = np.zeros((num_frames, num_rad))
for i, frame in enumerate(frame_range):
    density = (fromfile("gasdens%d.dat" % frame).reshape(num_rad, num_theta))
    averagedDensity = np.average(density, axis = 1) / surface_density
    averagedDensityCube[i] = averagedDensity

## Save Pickle File ##
storage = {}
storage['cube'] = averagedDensityCube
storage['scale_height'] = scale_height
storage['rad'] = rad
storage['id'] = fargo_par

random_id = random.randint(0, 9999)
pickle.dump(storage, open("averagedDensities%d.p" % random_id, "wb"))


