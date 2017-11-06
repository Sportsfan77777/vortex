"""
utility module for azimuthal-related plots
"""

import numpy as np
from scipy.signal import gaussian

import util

###############################################################################

### Helper Methods ###
def my_searchsorted(array, target):
    """ np.searchsorted, but it works all the time """
    for i, x in enumerate(array):
        if target > x:
            pass
        else:
            return i
    return len(array)

def get_radial_peak(averagedDensity, fargo_par):
    """ find peak in azimuthally-averaged density in the outer disk (i.e. the vortex) """
    ######## Get Parameters #########
    rad = fargo_par["rad"]

    ########### Method ##############
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    peak_rad_outer_index = np.argmax(averagedDensity[outer_disk_start : outer_disk_end])

    peak_index = outer_disk_start + peak_rad_outer_index
    peak_rad = rad[peak_index]
    peak_density = averagedDensity[peak_index]

    return peak_rad, peak_density

def get_radial_min(averagedDensity, peak_rad, fargo_par):
    """ find minimum in azimuthally-averaged density in the outer disk (i.e. the gap edge) """
    ######## Get Parameters #########
    rad = fargo_par["rad"]

    ########### Method ##############
    try:
        outer_disk_start = np.searchsorted(rad, 1.0) # look for min radial density beyond r = 1.1
        outer_disk_end = np.searchsorted(rad, peak_rad)
        min_rad_outer_index = np.argmin(averagedDensity[outer_disk_start : outer_disk_end])

        min_index = outer_disk_start + min_rad_outer_index
        min_rad = rad[min_index]
        min_density = averagedDensity[min_index]

        #print "Min", min_rad, min_density
        return min_rad, min_density
    except:
        # No Gap Yet
        return peak_rad, 0

def get_azimuthal_peak(density, fargo_par):
    """ return shift needed to shift vortex peak to 180 degrees """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    ########### Method ##############
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    argmax = np.argmax(density_segment)
    arg_r, arg_phi = np.unravel_index(argmax, np.shape(density_segment))

    ### Calculate shift for true center to 180 degrees ###
    middle = np.searchsorted(theta, np.pi)
    shift_peak = int(middle - arg_phi)

    return shift_peak

def get_azimuthal_center(density, fargo_par, threshold = 0.05):
    """ return shift needed to shift vortex center to 180 degrees """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############

    ### Identify center using threshold ###
    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, 1.1) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, 2.3) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    # Get peak in azimuthal profile
    avg_density = np.average(density_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_density)
    arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
    peak_rad = rad[arg_peak]

    # Zoom in on peak --- Average over half a scale height
    half_width = 0.25 * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    density_sliver = density[zoom_start : zoom_end]
    length = len(density_sliver); std = length / 3.0
    weights = gaussian(length, std)
    avg_density_sliver = np.average(density_sliver, weights = weights, axis = 0) # avg over rad

    # Move Minimum to Zero Degrees (vortex cannot cross zero)
    arg_min = np.argmin(avg_density_sliver)
    shift_min = int(0 - arg_min)
    avg_density_sliver = np.roll(avg_density_sliver, shift_min)

    # Spot two threshold crossovers
    left_edge = my_searchsorted(avg_density_sliver, threshold)
    right_edge = len(theta) - my_searchsorted(avg_density_sliver[::-1], threshold) - 1

    center = (left_edge + right_edge) / 2.0

    ### Calculate shift for true center to 180 degrees ###
    middle = np.searchsorted(theta, np.pi)
    shift_c = int(middle - (center - shift_min))

    return shift_c

### Data ###

def get_profiles(density, fargo_par, args, normalize = True, shift = None):
    """ Gather azimuthal radii and profiles """
    
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############

    if normalize:
        density /= surface_density_zero

    if shift is not None:
        density = np.roll(density, shift_c)

    # Find Peak in Radial Profile (in Outer Disk)
    averagedDensity = np.average(density, axis = 1)
    peak_rad, peak_density = get_radial_peak(averagedDensity, fargo_par)

    # Gather Azimuthal Profiles
    num_profiles = args.num_profiles
    spread = (args.num_scale_heights / 2.0) * scale_height

    azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
    azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
    azimuthal_profiles = [density[azimuthal_index, :] for azimuthal_index in azimuthal_indices]

    return azimuthal_radii, azimuthal_profiles

###############################################################################