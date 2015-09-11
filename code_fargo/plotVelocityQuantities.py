"""
plots velocityField as a vector field and as a magnitude field
plots vorticity
(both centered on planet)

ideally, plots in polar and cartesian
then, makes movies
"""

import math
import numpy as np

from matplotlib import rc
from matplotlib import pyplot as plot

from pylab import fromfile

import subprocess

# use LaTeX, choose nice some looking fonts and tweak some settings
"""
rc('font', family='serif')
rc('font', size=16)
rc('legend', fontsize=16)
rc('legend', numpoints=1)
rc('legend', handlelength=1.5)
rc('legend', frameon=False)
rc('xtick.major', pad=7)
rc('xtick.minor', pad=7)
rc('text', usetex=True)
rc('text.latex', preamble=[r'\usepackage[T1]{fontenc}',
                        r'\usepackage{amsmath}',
                        r'\usepackage{txfonts}',
                        r'\usepackage{textcomp}'])
"""

# Plot Parameters
cmap = "Oranges_r"

# Curl function
def curl(v_rad, v_theta, rad, theta):
    """ z-component of the curl (because this is a 2-D simulation)"""
    ### Start Differentials ###
    # d_r
    rad_shift = np.roll(rad, 1)
    dr = (rad - rad_shift)[1:]

    # d_t
    theta_shift = np.roll(theta, 1)
    dt = (theta - theta_shift)[1:]

    # dv_rad
    v_rad_shift = np.roll(v_rad, 1, axis = 1) # something is wrong here
    dv_rad = (v_rad - v_rad_shift)[:, 1:]

    # dv_theta
    v_theta_shift = np.roll(v_theta, 1, axis = 1) # something is wrong here
    dv_theta = (v_theta - v_theta_shift)[:, 1:]

    ### End Differentials ###

    # z-Determinant
    partial_one = rad[:-1, None] * dv_theta[:-1] / dr[:, None] # Note: dr is one shorter than rad
    partial_two = dv_rad[:-1, :] / dt[None, :] # Note: dt is d_theta, not d_time!!!!!!!!!

    # Source: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
    z_curl = (partial_one - partial_two) / rad[:-1, None]
    return z_curl

# Get Parameters
num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

for i in range(105):
    density = fromfile("gasdens%d.dat" % i).reshape(num_rad, num_theta)
    v_rad = fromfile("gasvrad%d.dat" % i).reshape(num_rad, num_theta)
    v_theta = fromfile("gasvtheta%d.dat" % i).reshape(num_rad, num_theta)

    # Roll Theta (so planet is not centered at theta = zero)
    roll = False
    if roll:
        theta = np.roll(theta, int(num_theta / 2))
        v_theta = np.roll(v_theta, int(num_theta / 2), axis = 1)

    # Velocity
    v_x = v_rad * np.cos(theta) - v_theta * np.sin(theta)
    v_y = v_rad * np.sin(theta) + v_theta * np.cos(theta)
    
    # Vorticity
    w = curl(v_rad, v_theta, rad, theta)
    
    ###### (1) Plot Velocity Field and Magnitude in Polar ######
    plot_velocity = True # change to option parser '-v'

    if plot_velocity:
        ### (1a.i) Field in Polar ###
        fig = plot.figure()
        ax = fig.add_subplot(111, polar = True)
        subset_r = range(int(num_rad / 2.7), int(num_rad / 1.5), 2)
        subset_t = range(0, int(num_theta), 10)
        ax.quiver(theta[subset_t], rad[subset_r], 
                  (v_x[subset_r])[:,subset_t], (v_y[subset_r])[:,subset_t])

        fontsize = 14

        plot.title("Velocity Field at Timestep %d" % i, fontsize = fontsize + 2)
        plot.savefig("velocityField/frame%03d_polar.png" % i, bbox_inches = 'tight')
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

        ### (1a.ii) Field in "r-theta Cartesian" ###
        fig = plot.figure()
        ax = fig.add_subplot(111)
        subset_r = range(int(num_rad / 2.7), int(num_rad / 1.5))
        subset_t = range(0, int(num_theta))
        ax.quiver(theta[subset_t], rad[subset_r], 
                  (v_rad[subset_r])[:,subset_t], (v_theta[subset_r])[:,subset_t])

        fontsize = 14

        plot.title("Velocity Field at Timestep %d" % i, fontsize = fontsize + 1)
        plot.savefig("velocityField/frame%03d_cart.png" % i, bbox_inches = 'tight')
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

        ### (1b.i) Magnitude in Polar ###
        fig = plot.figure()
        ax = fig.add_subplot(111, polar = True)
        subset_r = range(int(num_rad / 4.5), int(num_rad / 1.2), 2)
        subset_t = range(0, int(num_theta), 10)
        ax.pcolormesh(theta[subset_t], rad[subset_r], 
                     np.power((v_x[subset_r])[:,subset_t], 2) 
                     + np.power((v_y[subset_r])[:,subset_t], 2), cmap = cmap)

        fontsize = 14

        plot.title("Velocity Magnitude at Timestep %d" % i, fontsize = fontsize + 2)
        plot.savefig("velocityMagnitude/frame%03d_polar.png" % i, bbox_inches = 'tight')
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)

        ### (1b.ii) Magnitude in "r-theta Cartesian" ###
        fig = plot.figure()
        ax = fig.add_subplot(111)
        subset_r = range(int(num_rad / 4.5), int(num_rad / 1.2))
        subset_t = range(0, int(num_theta))
        ax.pcolormesh(theta[subset_t], rad[subset_r], 
                     np.power(v_rad[subset_r], 2) + np.power(v_theta[subset_r], 2), cmap = cmap)

        fontsize = 14

        plot.title("Velocity Magnitude at Timestep %d" % i, fontsize = fontsize + 1)
        plot.savefig("velocityMagnitude/frame%03d_cart.png" % i, bbox_inches = 'tight')
        #plot.show()
        plot.close(fig) # Close Figure (to avoid too many figures)


    ###### (2) Plot Velocity Field and Magnitude in Polar ######
    plot_vorticity = True # change to option parse '-w'

    if plot_vorticity:
        ### (2.i) Magnitude in Polar ###
        fig = plot.figure()
        ax = fig.add_subplot(111, polar = True)
        subset_r = range(int(num_rad / 4.5), int(num_rad / 1.2))
        subset_t = range(0, int(num_theta) - 1)
        ax.pcolormesh(theta[subset_t], rad[subset_r], (w[subset_r])[:,subset_t], cmap = cmap)

        fontsize = 14

        plot.title("Vorticity at Timestep %d" % i, fontsize = fontsize + 2)
        plot.savefig("vorticity/frame%03d_polar.png" % i, bbox_inches = 'tight')

        plot.close(fig)

        ### (2.ii) Magnitude in "r-theta Cartesian" ###
        fig = plot.figure()
        ax = fig.add_subplot(111)
        subset_r = range(int(num_rad / 4.5), int(num_rad / 1.2))
        subset_t = np.concatenate((range(0, int(num_theta / 8)), 
                            range(int(7 * num_theta / 8), int(num_theta - 1))))
        ax.pcolormesh(theta[subset_t], rad[subset_r], (w[subset_r])[:,subset_t], cmap = cmap)

        fontsize = 14

        plot.title("Vorticity at Timestep %d" % i, fontsize = fontsize + 2)
        plot.savefig("vorticity/frame%03d_cart.png" % i, bbox_inches = 'tight')

        plot.close(fig)

##### Make Movies #####

if plot_velocity:
    command = "avconv -framerate 5 -f image2 -vf scale=-2:720 -i velocityField/frame%03d_cart.png -b 65536k velocityField/velocityField_cart.mov"
    split_command = command.split()
    subprocess.Popen(split_command)

    command = "avconv -framerate 5 -f image2 -vf scale=-2:720 -i velocityMagnitude/frame%03d_cart.png -b 65536k velocityMagnitude/velocityMagnitude_cart.mov"
    split_command = command.split()
    subprocess.Popen(split_command)

    command = "avconv -framerate 5 -f image2 -i velocityField/frame%03d_polar.png -b 65536k velocityField/velocityField_polar.mov"
    split_command = command.split()
    subprocess.Popen(split_command)

    command = "avconv -framerate 5 -f image2 -i velocityMagnitude/frame%03d_polar.png -b 65536k velocityMagnitude/velocityMagnitude_polar.mov"
    split_command = command.split()
    subprocess.Popen(split_command)

if plot_vorticity:
    command = "avconv -framerate 5 -f image2 -i -vf scale=-2:720 vorticity/frame%03d_cart.png -b 65536k vorticity/vorticity_cart.mov"
    split_command = command.split()
    subprocess.Popen(split_command)

    command = "avconv -framerate 5 -f image2 -i vorticity/frame%03d_polar.png -b 65536k vorticity/vorticity_polar.mov"
    split_command = command.split()
    subprocess.Popen(split_command)
