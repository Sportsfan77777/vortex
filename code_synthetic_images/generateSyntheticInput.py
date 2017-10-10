"""
Usage:
python generateSyntheticInput.py

Does
(1) Generates secondary input: radial.dat, azimuthal.dat, temperature.dat, grain.dat
(2) Convert density distributions from code units to cgs
(3) Options to carve out inner cavity (default: yes) and scale dust sizes (default: no)
(4) Centers vortices (using cm, hcm, and mm peaks as a guide -- For T = 1000, center assumed to be fixed)
(5) Re-samples at a lower resolution
(6) Interpolates dust to continuous distribution of grain sizes (e.g. 100)
(7) Outputs "gasddens%d.dat", a text file input for the synthetic image generator (access with "ln -s gasddens%d.dat" density.dat)

Assumes directory structure with adjacent access to directories containing different grain sizes
"""

import sys
import os
import subprocess
import pickle
import glob
from multiprocessing import Pool

import math
import numpy as np
from scipy import interpolate as sp_int

from pylab import fromfile

directories = ["cm", "hcm", "mm", "hmm", "hum", "um"]





