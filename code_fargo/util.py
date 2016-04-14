"""
stores all function calls

Usage:
not to be called
"""

import os
import glob

import utilTorque

#### Miscellaneous ####

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


#### Torque Functions ####

""" """
torque = utilTorque.torque

inner_torque = utilTorque.inner_torque

outer_torque = utilTorque.outer_torque

inner_torque_contributions = utilTorque.inner_torque_contributions

outer_torque_contributions = utilTorque.outer_torque_contributions

#### Vorticity ####

