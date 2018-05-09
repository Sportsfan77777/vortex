"""

Usage:
not for direct calls
"""

import numpy as np

# Curl function
def velocity_curl(v_rad, v_theta, rad, theta, average = False, rossby = True, not_rotating_frame = False):
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

    if average:
        z_curl = (partial_one[:, 1:]) / rad[1:, None]
    else:
        z_curl = (partial_one[:, 1:] - partial_two[1:, :]) / rad[1:, None]

    # Shift out of rotating frame (http://arxiv.org/pdf/astro-ph/0605237v2.pdf)
    if not_rotating_frame:
        z_curl += 2

    # Convert to Rossby number (divide by angular velocity)
    if rossby:
        omega = np.pow(rad, -1.5)
        z_curl /= omega[1:, None]

    return z_curl