"""
stores all function calls

Usage:
not to be called
"""

import os
import subprocess
import glob
import time
import pickle

import numpy as np

import utilTorque
import utilVorticity

#### Miscellaneous ####

def pickle_parameters():
    param_fn = "params.p"
    if not os.path.exists(param_fn):
        command = "python pickleParameters.py"
        split_command = command.split()
        subprocess.Popen(split_command)
    else:
        return False
    # Give it time to execute
    time.sleep(1) # one sec
    return True

def find_max_frame():
    # glob the right set of files
    fargo_fn = "fargo2D1D"
    if os.path.exists(fargo_fn):
        density_files = glob.glob("gasdens1D*.dat")
        number_index = 9 # the number starts at index 9 after 'gasdens1D'
    else:
        density_files = glob.glob("gasdens*.dat")
        number_index = 7 # the number starts at index 7 after 'gasdens'
    # Find the highest numbered file
    max_frame = 0
    for d_f in density_files:
        name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
        frame_number = int(name[number_index:]) # just 999
        if frame_number > max_frame:
            max_frame = frame_number
    return max_frame

def getWaveLocation(density, radii, theta, start_radius = 1.10, end_radius = 2.30):
    # Find Limits
    start_i = np.searchsorted(radii, start_radius)
    end_i = np.searchsorted(radii, end_radius)

    radii = radii[start_i : end_i] # truncate to selected range

    # Find maxes
    density_near_vortex = density[start_i : end_i]
    wave_locations = np.zeros(len(density_near_vortex))

    guess_theta = 2 * np.pi
    delta_theta = 10.0 * (np.pi / 180.0)
    for i, r_i in enumerate(radii):
        # Guess Bounds
        upper_guess = guess_theta # + delta_theta
        lower_guess = guess_theta - delta_theta # search only below previous location

        # Mask Values Away From Wave 
        density_ring = density_near_vortex[i]
        mask = np.ones(len(density_ring)) # 1 is masked

        if lower_guess < 0.0:
            lower_guess_i = np.searchsorted(theta, lower_guess + 2 * np.pi)
            upper_guess_i = np.searchsorted(theta, upper_guess)

            mask[lower_guess_i : ] = 0
            mask[ : upper_guess_i] = 0

        elif upper_guess > 2 * np.pi:
            lower_guess_i = np.searchsorted(theta, lower_guess)
            upper_guess_i = np.searchsorted(theta, upper_guess - 2 * np.pi)

            mask[lower_guess_i : ] = 0
            mask[ : upper_guess_i] = 0

        else:
            lower_guess_i = np.searchsorted(theta, lower_guess)
            upper_guess_i = np.searchsorted(theta, upper_guess)
            
            mask[lower_guess_i : upper_guess_i] = 0

        masked_density_ring = np.ma.array(density_ring, mask = mask)
        wave_index = np.ma.argmax(masked_density_ring, fill_value = -1.0)

        # Store Result + Update Next Guess
        wave_locations[i] = theta[wave_index]
        guess_theta = wave_locations[i]

    return radii, wave_locations

#### Torque Functions ####

""" """
torque = utilTorque.torque

inner_torque = utilTorque.inner_torque

outer_torque = utilTorque.outer_torque

inner_torque_contributions = utilTorque.inner_torque_contributions

outer_torque_contributions = utilTorque.outer_torque_contributions

#### Vorticity ####

velocity_curl = utilVorticity.velocity_curl

