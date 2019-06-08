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

def get_polar_grid(xs, ys, min_r = 0):
    """ return polar grid associated with cartesian to polar"""
    max_r = np.sqrt(np.power(np.max(xs), 2) + np.power(np.max(ys), 2))
    num_r = len(xs); num_theta = len(ys)

    # Set up rt-grid
    rs = np.linspace(min_r, max_r, num_r)
    thetas = np.linspace(0, 2 * np.pi, num_theta)

    rs_grid, thetas_grid = np.meshgrid(rs, thetas)
    return rs, thetas, rs_grid, thetas_grid

def polar_to_cartesian(data, rs, thetas, order = 3):
    """ Input: Data, Radial array for Data, Azimuthal array for Data"""
    # Source: http://stackoverflow.com/questions/2164570/reprojecting-polar-to-cartesian-grid
    # Overall, use an interpolator to space a proper grid instead of just converting the (r,t) coordinates to (x,y)

    # Note: Reversing thetas is necessary to maintain counter-clockwise rotation
    thetas = thetas[::-1]

    # Set up xy-grid
    xs, ys, xs_grid, ys_grid = get_cartesian_grid(rs)

    # Interpolators for the rt-grid (Functions!!!!)
    rs_interpolator = interpolate(rs, np.arange(len(rs)), bounds_error = False)
    thetas_interpolator = interpolate(thetas, np.arange(len(thetas)), bounds_error = False)

    # Match up xy-grid with rt-grid
    new_rs = np.sqrt(np.power(xs_grid, 2) + np.power(ys_grid, 2))
    new_thetas = np.arctan2(ys_grid, xs_grid)

    # Convert from [-pi, pi] to [0, 2 * pi]
    new_thetas[new_thetas < 0] = new_thetas[new_thetas < 0] + 2 * np.pi

    # Get new_rs and new_thetas with the interpolators
    new_interpolated_rs = rs_interpolator(new_rs.ravel())
    new_interpolated_thetas = thetas_interpolator(new_thetas.ravel())

    # Fix Bounds (outside of polar image, but inside cartesian image)
    new_interpolated_rs[new_rs.ravel() > max(rs)] = len(rs) - 1
    new_interpolated_rs[new_rs.ravel() < min(rs)] = 0

    cart_data = map_coordinates(data, np.array([new_interpolated_rs, new_interpolated_thetas]), order = order).reshape(new_rs.shape)

    return xs, ys, xs_grid, ys_grid, cart_data

def cartesian_to_polar(data, xs, ys, order = 3):
    # Set up rt-grid
    rs, thetas, rs_grid, thetas_grid = get_polar_grid(xs, ys)

    # Interpolators for the rt-grid (Functions!!!!)
    xs_interpolator = interpolate(xs, np.arange(len(xs)), bounds_error = False)
    ys_interpolator = interpolate(ys, np.arange(len(ys)), bounds_error = False)

    # Match up rt-grid with xy-grid
    new_xs = rs_grid * np.cos(thetas_grid)
    new_ys = rs_grid * np.sin(thetas_grid)

    # Get new_xs and new_ys with the interpolators
    new_interpolated_xs = xs_interpolator(new_xs.ravel())
    new_interpolated_ys = xs_interpolator(new_ys.ravel())

    polar_data = map_coordinates(data, np.array([new_interpolated_xs, new_interpolated_ys]), order = order).reshape(new_xs.shape)

    return rs, thetas, rs_grid, thetas_grid, polar_data

