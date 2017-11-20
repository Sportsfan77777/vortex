"""
stores utility functions for intensity

Usage:
not to be called
"""

import os
import pickle

import numpy as np


def record_contrast(intensity, xs, ys):
	""" Change this to assume polar coordinates """
    # Get intensity along x = 0
    zero_i = np.searchsorted(xs, 0)
    two_sliver = intensity[zero_i - 1 : zero_i + 1, :] # two columns around x = 0
    sliver = np.average(two_sliver, axis = 0)

    # Mask inner disk
    lower_i = np.searchsorted(ys, -1)
    upper_i = np.searchsorted(ys, 1)
    sliver[lower_i : upper_i] = 0.0

    # Find argmax (y-coor, and opposite y-coor)
    max_yi = np.argmax(sliver)
    max_y = ys[max_yi]

    opposite_i = np.searchsorted(ys, -max_y)

    # Mark contrast, max, opposite
    maximum = np.max(sliver)
    opposite = sliver[opposite_i]
    contrast = maximum / opposite

    return contrast, maximum, opposite

