"""
utility module for azimuthal-related plots
"""

import pickle

import numpy as np
from scipy.signal import gaussian

import util

###############################################################################

outer_start = 1.1
outer_end = 2.3

### Helper Methods ###

def my_searchsorted(array, target):
    """ np.searchsorted, but it works all the time """
    for i, x in enumerate(array):
        if target > x:
            pass
        else:
            return i
    return len(array)

def get_radial_peak(averagedDensity, fargo_par, start = outer_start, end = outer_end):
    """ find peak in azimuthally-averaged density in the outer disk (i.e. the vortex) """
    ######## Get Parameters #########
    rad = fargo_par["rad"]

    ########### Method ##############
    outer_disk_start = np.searchsorted(rad, start) # look for max radial density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
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

def get_peak(density, fargo_par, start = outer_start, end = outer_end):
    """ return location of peak in data in outer disk """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    ########### Method ##############
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
    density_segment = density[outer_disk_start : outer_disk_end]

    argmax = np.argmax(density_segment)
    arg_r, arg_phi = np.unravel_index(argmax, np.shape(density_segment))

    arg_r += outer_disk_start

    return arg_r, arg_phi

def get_azimuthal_peak(density, fargo_par, start = outer_start, end = outer_end):
    """ return shift needed to shift vortex peak to 180 degrees """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    ########### Method ##############
    arg_r, arg_phi = get_peak(density, fargo_par, start = start, end = end)

    ### Calculate shift for true center to 180 degrees ###
    middle = np.searchsorted(theta, np.pi)
    shift_peak = int(middle - arg_phi)

    return shift_peak

def get_azimuthal_center(density, fargo_par, threshold = 0.05, start = outer_start, end = outer_end):
    """ return shift needed to shift vortex center to 180 degrees """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############

    ### Identify center using threshold ###
    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
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

def shift_away_from_minimum(density, fargo_par, start = outer_start, end = outer_end):
    """ return shift needed to shift vortex center to 180 degrees """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############

    ### Identify center using threshold ###
    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
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

    return shift_min

def find_shift(density, reference_density, fargo_par, center = True, num_scale_heights = 4.0, min_shift = 0, max_shift = 359, num_shifts = 360, start = outer_start, end = outer_end):
    """ return shift needed to overlap the density with the reference_density """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############

    ### Get Slivers ###

    def get_sliver(density):
        ### Identify center using threshold ###
        # Search outer disk only
        outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
        outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
        density_segment = density[outer_disk_start : outer_disk_end]

        # Get peak in azimuthal profile
        avg_density = np.average(density_segment, axis = 1) # avg over theta
        segment_arg_peak = np.argmax(avg_density)
        arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
        peak_rad = rad[arg_peak]

        # Zoom in on peak
        half_width = 0.5 * num_scale_heights * scale_height
        zoom_start = np.searchsorted(rad, peak_rad - half_width)
        zoom_end = np.searchsorted(rad, peak_rad + half_width)

        density_sliver = density[zoom_start : zoom_end]
        return density_sliver

    density_sliver = get_sliver(density)
    reference_density_sliver = get_sliver(reference_density)

    ### Shift reference density away from the min
    if center:
        shift_away = shift_away_from_minimum(reference_density, fargo_par)
        reference_density_sliver = np.roll(reference_density_sliver, shift_away)

    ### Test shifts ###
    theta_degrees = theta * (180.0 / np.pi)

    mass_differences = np.zeros(num_shifts)
    possible_shifts = np.linspace(min_shift, max_shift, num_shifts)

    for i, possible_shift_i in enumerate(possible_shifts):
        shift = np.searchsorted(theta_degrees, possible_shift_i)
        density_sliver_tmp = np.roll(density_sliver, shift)

        diff = np.abs(reference_density_sliver - density_sliver_tmp)
        mass_differences[i] = np.sum(diff)

    ### The correct shift has the lowest mass difference ###
    shift_i = np.argmin(mass_differences)
    theta_shift = possible_shifts[shift_i]

    shift_to_roll_i = np.searchsorted(theta_degrees, theta_shift) # this is an index

    return shift_to_roll_i, theta_shift

def get_lookup_shift(frame, directory = "."):
    """ lookup the shift determined with findShifts.py """
    frames = pickle.load(open('%s/frame_lookup.p' % directory, 'r'))
    shifts = pickle.load(open('%s/shift_lookup.p' % directory, 'r'))

    frame_i = np.searchsorted(frames, frame) # Note: not an int!
    return int(round(shifts[frame_i], 0))

### Extract Values ###

def get_polar_contrast(data, fargo_par):
    """ for polar data, returns contrast between peak and opposite point """
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    ########### Method ##############
    # Get indices
    arg_r, arg_phi = get_peak(data, fargo_par)

    # Get index of opposite
    phi = theta[arg_phi]
    opposite_phi = (phi + np.pi) % (2 * np.pi)
    arg_opposite = np.searchsorted(theta, opposite_phi)

    # Get target values
    data_peak = data[arg_r, arg_phi]
    data_opposite = data[arg_r, arg_opposite]

    contrast = data_peak / data_opposite

    return contrast, data_peak, data_opposite

def get_contrast(data, fargo_par, normalize = False, threshold = 0.5, sliver_width = 0.5, start = outer_start, end = outer_end):
    """ Get azimuthal extent at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]

    ########### Method ##############

    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
    data_segment = data[outer_disk_start : outer_disk_end]

    # Get peak in azimuthal profile
    avg_data = np.average(data_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_data)
    arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
    peak_rad = rad[arg_peak]

    # Zoom in on peak --- Average over half a scale height
    half_width = (0.5 * sliver_width) * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    data_sliver = data[zoom_start : zoom_end]
    length = len(data_sliver); std = length / 3.0
    weights = gaussian(length, std)
    azimuthal_profile = np.average(data_sliver, weights = weights, axis = 0) # avg over rad to get azimuthal profile

    if normalize:
        azimuthal_profile /= np.max(azimuthal_profile)

    contrast = np.max(azimuthal_profile) / np.min(azimuthal_profile)
    return contrast

def get_azimuthal_profile(data, fargo_par, normalize = False, threshold = 0.5, sliver_width = 0.5, start = outer_start, end = outer_end):
    """ Get azimuthal extent at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]

    ########### Method ##############

    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
    data_segment = data[outer_disk_start : outer_disk_end]

    # Get peak in azimuthal profile
    avg_data = np.average(data_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_data)
    arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
    peak_rad = rad[arg_peak]

    # Zoom in on peak --- Average over half a scale height
    half_width = (0.5 * sliver_width) * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    data_sliver = data[zoom_start : zoom_end]
    length = len(data_sliver); std = length / 3.0
    weights = gaussian(length, std)
    azimuthal_profile = np.average(data_sliver, weights = weights, axis = 0) # avg over rad to get azimuthal profile

    return azimuthal_profile

def get_mean_azimuthal_profile(data, fargo_par, normalize = False, threshold = 0.5, sliver_width = 1.5, start = outer_start, end = outer_end):
    """ Get azimuthal extent at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]

    ########### Method ##############

    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
    data_segment = data[outer_disk_start : outer_disk_end]

    # Get peak in azimuthal profile
    avg_data = np.average(data_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_data)
    arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
    peak_rad = rad[arg_peak]

    # Zoom in on peak --- Average over half a scale height
    half_width = (0.5 * sliver_width) * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    data_sliver = data[zoom_start : zoom_end] # This one is the mean! The other one is Gaussian-weighted!
    mean_azimuthal_profile = np.average(data_sliver, axis = 0) # avg over rad to get azimuthal profile

    return mean_azimuthal_profile

def get_extent(data, fargo_par, normalize = False, threshold = 0.5, sliver_width = 0.5, start = outer_start, end = outer_end):
    """ Get azimuthal extent at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]

    ########### Method ##############

    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
    data_segment = data[outer_disk_start : outer_disk_end]

    # Get peak in azimuthal profile
    avg_data = np.average(data_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_data)
    arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
    peak_rad = rad[arg_peak]

    # Zoom in on peak --- Average over half a scale height
    half_width = (0.5 * sliver_width) * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    data_sliver = data[zoom_start : zoom_end]
    length = len(data_sliver); std = length / 3.0
    weights = gaussian(length, std)
    azimuthal_profile = np.average(data_sliver, weights = weights, axis = 0) # avg over rad to get azimuthal profile

    if normalize:
        azimuthal_profile /= np.max(azimuthal_profile)

    # Move minimum to theta = zero
    arg_min = np.argmin(azimuthal_profile)
    shift_min = int(0 - arg_min)
    azimuthal_profile = np.roll(azimuthal_profile, shift_min)

    # Find extents with the threshold
    left_theta_i = my_searchsorted(azimuthal_profile, threshold)
    right_theta_i = len(theta) - (my_searchsorted(azimuthal_profile[::-1], threshold)) - 1

    left_theta = theta[left_theta_i]
    right_theta = theta[right_theta_i]

    extent = right_theta - left_theta
    return extent

def get_azimuthal_bounds(data, fargo_par, radial_center = None, normalize = False, threshold = 0.5, sliver_width = 0.5, start = outer_start, end = outer_end):
    """ Get azimuthal edges at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]

    ########### Method ##############

    if radial_center is None:
        # Search outer disk only
        outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
        outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
        data_segment = data[outer_disk_start : outer_disk_end]

        # Get peak in azimuthal profile
        avg_data = np.average(data_segment, axis = 1) # avg over theta
        segment_arg_peak = np.argmax(avg_data)
        arg_peak = np.searchsorted(rad, rad[outer_disk_start + segment_arg_peak])
        peak_rad = rad[arg_peak]
    else:
        peak_rad = radial_center

    # Zoom in on peak --- Average over half a scale height
    half_width = (0.5 * sliver_width) * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    data_sliver = data[zoom_start : zoom_end]
    length = len(data_sliver); std = length / 3.0
    weights = gaussian(length, std)
    azimuthal_profile = np.average(data_sliver, weights = weights, axis = 0) # avg over rad to get azimuthal profile

    if normalize:
        azimuthal_profile /= np.max(azimuthal_profile)

    # Move minimum to theta = zero
    arg_min = np.argmin(azimuthal_profile)
    shift_min = int(0 - arg_min)
    azimuthal_profile = np.roll(azimuthal_profile, shift_min)

    # Find extents with the threshold
    left_theta_i = my_searchsorted(azimuthal_profile, threshold)
    right_theta_i = len(theta) - (my_searchsorted(azimuthal_profile[::-1], threshold)) - 1

    left_theta = theta[left_theta_i]
    right_theta = theta[right_theta_i]

    return left_theta, right_theta

def get_radial_extent(data, fargo_par, normalize = False, threshold = 0.5, sliver_width = 20.0, start = 1.1, end = 2.5):
    """ Get radial extent at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]

    ########### Method ##############

    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
    rad_segment = rad[outer_disk_start : outer_disk_end]
    data_segment = data[outer_disk_start : outer_disk_end]

    # Move minimum to theta = zero (first get peak in azimuthally-averaged profile)
    avg_data = np.average(data_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_data)
    arg_peak = np.searchsorted(rad_segment, rad_segment[segment_arg_peak])
    peak_rad = rad[arg_peak]

    arg_min = np.argmin(data_segment[arg_peak])
    shift_min = int(0 - arg_min)
    data_segment = np.roll(data_segment, shift_min, axis = -1)

    # Get peak in azimuthal profile
    avg_data = np.average(data_segment, axis = 0) # avg over rad
    arg_peak = np.argmax(avg_data)
    peak_theta = theta[arg_peak]

    # Zoom in on center --- Average over sliver width
    half_width = (0.5 * sliver_width) * (np.pi / 180.0)
    zoom_start = np.searchsorted(theta, peak_theta - half_width)
    zoom_end = np.searchsorted(theta, peak_theta + half_width)

    data_sliver = data_segment[:, zoom_start : zoom_end]
    length = len(data_sliver[0]); std = length / 3.0
    weights = gaussian(length, std)
    radial_profile = np.average(data_sliver, weights = weights, axis = 1) # avg over theta to get radial profile

    if normalize:
        radial_profile /= np.max(radial_profile)

    # Find extents with the threshold
    left_rad_i = my_searchsorted(radial_profile, threshold)
    right_rad_i = len(rad_segment) - (my_searchsorted(radial_profile[::-1], threshold)) - 1

    left_rad = rad_segment[left_rad_i]
    right_rad = rad_segment[right_rad_i]

    extent = right_rad - left_rad
    middle = left_rad + 0.5 * (extent)
    return extent, middle

def get_radial_bounds(data, fargo_par, normalize = False, threshold = 0.5, sliver_width = 20.0, start = 1.1, end = 2.5):
    """ Get radial extent at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]

    ########### Method ##############

    # Search outer disk only
    outer_disk_start = np.searchsorted(rad, start) # look for max density beyond r = 1.1
    outer_disk_end = np.searchsorted(rad, end) # look for max density before r = 2.3
    rad_segment = rad[outer_disk_start : outer_disk_end]
    data_segment = data[outer_disk_start : outer_disk_end]

    # Move minimum to theta = zero (first get peak in azimuthally-averaged profile)
    avg_data = np.average(data_segment, axis = 1) # avg over theta
    segment_arg_peak = np.argmax(avg_data)
    arg_peak = np.searchsorted(rad_segment, rad_segment[segment_arg_peak])
    peak_rad = rad[arg_peak]

    arg_min = np.argmin(data_segment[arg_peak])
    shift_min = int(0 - arg_min)
    data_segment = np.roll(data_segment, shift_min, axis = -1)

    # Get peak in azimuthal profile
    avg_data = np.average(data_segment, axis = 0) # avg over rad
    arg_peak = np.argmax(avg_data)
    peak_theta = theta[arg_peak]

    # Zoom in on center --- Average over sliver width
    half_width = (0.5 * sliver_width) * (np.pi / 180.0)
    zoom_start = np.searchsorted(theta, peak_theta - half_width)
    zoom_end = np.searchsorted(theta, peak_theta + half_width)

    data_sliver = data_segment[:, zoom_start : zoom_end]
    length = len(data_sliver[0]); std = length / 3.0
    weights = gaussian(length, std)
    radial_profile = np.average(data_sliver, weights = weights, axis = 1) # avg over theta to get radial profile

    if normalize:
        radial_profile /= np.max(radial_profile)

    # Find extents with the threshold
    left_rad_i = my_searchsorted(radial_profile, threshold)
    right_rad_i = len(rad_segment) - (my_searchsorted(radial_profile[::-1], threshold)) - 1

    left_rad = rad_segment[left_rad_i]
    right_rad = rad_segment[right_rad_i]

    return left_rad, right_rad

### Data ###

def get_profiles(density, fargo_par, args, normalize = False, shift = None, start = outer_start, end = outer_end):
    """ Gather azimuthal radii and profiles (doesn't have to be density) """
    
    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    scale_height = fargo_par["AspectRatio"]
    surface_density_zero = fargo_par["Sigma0"]

    ########### Method ##############

    if normalize:
        density /= surface_density_zero

    if shift is not None:
        density = np.roll(density, shift)

    # Find Peak in Radial Profile (in Outer Disk)
    averagedDensity = np.average(density, axis = 1)
    peak_rad, peak_density = get_radial_peak(averagedDensity, fargo_par, start = start, end = end)

    # Gather Azimuthal Profiles
    num_profiles = args.num_profiles

    if num_profiles == 1:
        azimuthal_index = np.searchsorted(rad, peak_rad)
        azimuthal_profile = density[azimuthal_index]
        
        return peak_rad, azimuthal_profile
    else:
        spread = (args.num_scale_heights / 2.0) * scale_height
        
        azimuthal_radii = np.linspace(peak_rad - spread, peak_rad + spread, num_profiles)
        azimuthal_indices = [np.searchsorted(rad, this_radius) for this_radius in azimuthal_radii]
        azimuthal_profiles = [density[azimuthal_index, :] for azimuthal_index in azimuthal_indices]

        return azimuthal_radii, azimuthal_profiles

###############################################################################

### Plotting ###

def get_max_y(size, taper_time):
    """ return size name corresponding to size number """
    max_ys = {}
    if taper_time < 10.1:
        max_ys[1.0] = 1500; max_ys[0.3] = 400; max_ys[0.1] = 100
        max_ys[0.03] = 40; max_ys[0.01] = 20; max_ys[0.0001] = 4.0

        return max_ys[size]
    else:
        max_ys[1.0] = 600; max_ys[0.3] = 125; max_ys[0.1] = 30
        max_ys[0.03] = 20; max_ys[0.01] = 10; max_ys[0.0001] = 2.5

        return max_ys[size]

###############################################################################

### Analytic ###

def get_analytic_profile(angle, r, dr, dtheta, aspect_ratio, S, max_density = 1, scale_height = 0.06):
    """ Calculates analytic azimuthal dust density profile for a given aspect ratio and S = St / \delta """

    def semiminor_axis(angle, dr, dtheta):
        """ maps angle in vortex to semiminor axis input to dust distribution (Eq. 65 from Lyra + Lin 13) """
        # Note: The azimuthal edge of the vortex should map to the radial half-width (in units of the scale height)
        return (dr / scale_height) * (angle / dtheta)

    #def semiminor_axis2(angle, dr, dtheta):
    #    """ old method: incorrect probably??? """
    #    r_over_dr = r / dr
    #    return r_over_dr * (angle * (np.pi / 180.0)) * 0.06 * 2

    def scale_function_sq(aspect_ratio):
        xi = 1 + aspect_ratio**(-2); vorticity = 1.5 / (aspect_ratio - 1)

        first_term = 2.0 * vorticity * aspect_ratio
        second_term = xi**(-1) * (2 * (vorticity)**2 + 3)

        return first_term - second_term

    x = semiminor_axis(angle, dr, dtheta)
    f_sq = scale_function_sq(aspect_ratio)

    coeff = max_density * (S + 1)**(1.5)
    exp = np.exp(-(S + 1) * x**2 * f_sq / 2.0)

    return coeff * exp