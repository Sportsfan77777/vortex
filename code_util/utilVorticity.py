"""

Usage:
not for direct calls
"""

import numpy as np

# Curl function
def velocity_curl(v_rad, v_theta, rad, theta, average = False, rossby = True, residual = True):
    """ z-component of the curl (because this is a 2-D simulation)"""
    # subtract rotating frame and add real Keplerian velocity
    keplerian_velocity = rad * (np.power(rad, -1.5) - 1) # in rotating frame, v_k = r * (r^-1.5 - r_p^-1.5)
    v_theta -= keplerian_velocity[:, None]

    v_keplerian_real = np.power(rad, -0.5)
    v_theta += v_keplerian_real[:, None]

    v_theta_average = np.average(v_theta, axis = 1)

    if residual:
        # subtract out the real Keplerian velocity
        #v_keplerian = np.power(rad, -0.5)
        v_theta -= v_keplerian_real[:, None]
        
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
    #z_curl += 2

    # Convert to Rossby number (divide by [2 * angular velocity])
    if rossby:
        omega = np.power(rad, -1.5) # Change this to v_theta / r ????????
        #omega = v_theta_average / rad
        z_curl /= (2.0 * omega[1:, None])

    return z_curl

def velocity_curl_3D(v_rad, v_theta, rad, theta, phi = 0, average = False, rossby = True, residual = True):
    """ z-component of the curl (because this is a 2-D simulation)"""
    phi = int(round(phi, 0))

    # subtract rotating frame and add real Keplerian velocity
    keplerian_velocity = rad * (np.power(rad, -1.5) - 1) # in rotating frame, v_k = r * (r^-1.5 - r_p^-1.5)
    v_theta -= keplerian_velocity[:, None]

    v_keplerian_real = np.power(rad, -0.5)
    v_theta += v_keplerian_real[:, None]

    v_theta_average = np.average(v_theta, axis = 1)

    if residual:
        # subtract out the real Keplerian velocity
        #v_keplerian = np.power(rad, -0.5)
        v_theta -= v_keplerian_real[:, None]
        
    ### Start Differentials ###
    d_rad = np.diff(rad)
    d_theta = theta[1] - theta[0]

    dv_rad = np.diff(v_rad, axis = 1)
    dv_theta = np.diff(rad[:, None] * v_theta, axis = 0)

    dv_rad = dv_rad[:, :, phi]
    dv_theta = dv_theta[:, :, phi]
    ### End Differentials ###

    # z-Determinant
    partial_one = dv_theta[:, 1:] / d_rad[None, :]
    partial_two = dv_rad / d_theta

    if average:
        z_curl = (partial_one[:, :]) / rad[None, 1:]
    else:
        z_curl = (partial_one[:, :] - partial_two[1:, :]) / rad[None, 1:]

    # Shift out of rotating frame (http://arxiv.org/pdf/astro-ph/0605237v2.pdf)
    #z_curl += 2

    # Convert to Rossby number (divide by [2 * angular velocity])
    if rossby:
        omega = np.power(rad, -1.5) # Change this to v_theta / r ????????
        #omega = v_theta_average / rad
        z_curl /= (2.0 * omega[None, 1:])

    return z_curl