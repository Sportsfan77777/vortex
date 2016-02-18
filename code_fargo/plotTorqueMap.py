"""
plots contributions to torque for each cell

Usage:
python plotTorqueMap.py
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

# Constants #
BigG = 1

def curl(r1, theta1, r2, theta2):
    pass


def torque(radius, theta, density):
	pass