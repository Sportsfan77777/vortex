"""
functions related to converting between cartesian and polar

Usage:
not to be called
"""

import numpy as np

from scipy import signal
from scipy.interpolate import interp1d as interpolate
from scipy.ndimage import map_coordinates

def get_cartesian_grid(rad):
    """ return cartesian grid associated with polar to cartesian"""
    max_r = rad[-1]
    resolution = 2 * len(rad)

    # Set up xy-grid
    xs = np.linspace(-max_r, max_r, resolution)
    ys = np.linspace(-max_r, max_r, resolution)

    xs_grid, ys_grid = np.meshgrid(xs, ys)
    return xs, ys, xs_grid, ys_grid


def polar_to_cartesian(data, rs, thetas, order = 3):
    """ Input: Data, Radial array for Data, Azimuthal array for Data"""
    # Source: http://stackoverflow.com/questions/2164570/reprojecting-polar-to-cartesian-grid

    # Note: Reversing thetas is necessary to maintain counter-clockwise rotation
    thetas = thetas[::-1]

    # Set up xy-grid
    xs, ys, xs_grid, ys_grid = get_cartesian_grid(rs)

    # Interpolate rt-grid

    interpolated_rs = interpolate(rs, np.arange(len(rs)), bounds_error = False)
    interpolated_thetas = interpolate(thetas, np.arange(len(thetas)), bounds_error = False)

    # Match up xy-grid with rt-grid

    new_rs = np.sqrt(np.power(xs_grid, 2) + np.power(ys_grid, 2))
    new_thetas = np.arctan2(ys_grid, xs_grid)

    # Convert from [-pi, pi] to [0, 2 * pi]
    new_thetas[new_thetas < 0] = new_thetas[new_thetas < 0] + 2 * np.pi

    new_interpolated_rs = interpolated_rs(new_rs.ravel())
    new_interpolated_thetas = interpolated_thetas(new_thetas.ravel())

    # Fix Bounds (outside of polar image, but inside cartesian image)
    new_interpolated_rs[new_rs.ravel() > max(rs)] = len(rs) - 1
    new_interpolated_rs[new_rs.ravel() < min(rs)] = 0

    cart_data = map_coordinates(data, np.array([new_interpolated_rs, new_interpolated_thetas]), order = order).reshape(new_rs.shape)

    return xs, ys, xs_grid, ys_grid, cart_data

