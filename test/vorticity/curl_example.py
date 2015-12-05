"""
test curl function
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


# Curl function
def curl(v_rad, v_theta, rad, theta):
    """ z-component of the curl (because this is a 2-D simulation)"""
    ### Start Differentials ###
    d_rad = np.diff(rad)
    d_theta = np.diff(theta)

    dv_rad = np.diff(v_rad, axis = 1)
    dv_theta = np.diff(rad[:, None] * v_theta, axis = 0)
    ### End Differentials ###

    # z-Determinant
    partial_one = dv_theta / d_rad[:, None]
    partial_two = dv_rad / d_theta

    z_curl = (partial_one[:, 1:] - partial_two[1:, :]) / rad[1:, None]
    return z_curl

#### DATA ####
num_rad = 3000
num_theta = 3000
rs = np.linspace(0.2, 4, num_rad)
thetas = np.linspace(0, 2 * np.pi, num_theta)

mesh = np.meshgrid(thetas, rs)
#print (mesh[1])[0,0], (mesh[0])[0,0]

#print np.shape(mesh)

def f_rad(r, theta):
	return 3 * r + np.cos(theta)

def f_theta(r, theta):
	return 2 * r**2 * (np.sin(theta))**2

def solution(r, theta):
	term_one = 6 * r * (np.sin(theta))**2
	term_two = np.sin(theta) / r
	return term_one + term_two

curl_solution = solution(mesh[1], mesh[0])

frad_mesh = f_rad(mesh[1], mesh[0])
ftheta_mesh = f_theta(mesh[1], mesh[0])

curl_calculation = curl(frad_mesh, ftheta_mesh, rs, thetas)

diff = curl_solution[1:, 1:] - curl_calculation

#### PLOTTING ####

fig = plot.figure()
ax = fig.add_subplot(111, polar = True)

result = ax.pcolormesh(thetas[1:], rs[1:], diff)
fig.colorbar(result)
result.set_clim(-0.01, 0.01)
#plot.savefig("notsaving.png")
plot.show()

def analyze(f):
	print np.median(f), np.min(f), np.max(f), f

print "Calculated"
analyze(curl_calculation)
print "Actual"
analyze(curl_solution)
print "Diff"
analyze(diff)

#print np.median(diff), np.max(diff)