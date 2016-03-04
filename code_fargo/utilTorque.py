"""

Usage:
not for direct calls
"""

import numpy as np

# Constants
BigG = 1

# Torque Function Element by Element
def torque(radius, theta, density, planet_mass = 0.005):
    """ 'torque' = 'r cross F' """
    # Differentials
    d_rad = np.diff(radius)
    d_theta = np.diff(theta)

    d_rad = np.append(d_rad, d_rad[-1])
    d_theta = np.append(d_theta, d_theta[-1])

    # Relevant Vectors + Quantities
    r_element = np.array([np.outer(radius, np.cos(theta)), np.outer(radius, np.sin(theta))]) # to fluid element 
    r_diff = np.array([np.outer(radius, np.cos(theta)) - 1, np.outer(radius, np.sin(theta))]) # to fluid element 

    dist_sq = np.einsum("ijk,ijk->jk", r_diff, r_diff)

    # Torque
    coeff = BigG * planet_mass * density / dist_sq
    direction = np.cross(r_element, r_diff, axis = 0)

    torque_density = coeff * direction
    area = radius[:, None] * np.outer(d_rad, d_theta)

    print area[-1, -1]
    print np.shape(area)

    return torque_density * area

def total_inner_torque(radius, theta, density, planet_loc = 1.0):
    """ sums up contribution from inner disk only """
    planet_index = np.searchsorted(radius, planet_loc)

    inner_torque = np.sum(density[:planet_index])
    return inner_torque

def total_outer_torque():
    """ sums up contribution from outer disk only """
    planet_index = np.searchsorted(radius, planet_loc)

    outer_torque = np.sum(density[planet_index:])
    return outer_torque
