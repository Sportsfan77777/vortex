"""
utility module for azimuthal-related plots
"""

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

### Extract Values ###

def get_contrast(data, fargo_par):
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


def get_extent(data, fargo_par, threshold = 0.5, sliver_width = 0.5, start = outer_start, end = outer_end):
    """ Get azimuthal extent at peak across a given threshold """

    ######## Get Parameters #########
    rad = fargo_par["rad"]
    theta = fargo_par["theta"]

    ########### Method ##############

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
    half_width = (0.5 * sliver_width) * scale_height
    zoom_start = np.searchsorted(rad, peak_rad - half_width)
    zoom_end = np.searchsorted(rad, peak_rad + half_width)

    density_sliver = density[zoom_start : zoom_end]
    length = len(density_sliver); std = length / 3.0
    weights = gaussian(length, std)
    azimuthal_profile = np.average(density_sliver, weights = weights, axis = 0) # avg over rad to get azimuthal profile

    # Move minimum to theta = zero
    arg_min = np.argmin(azimuthal_profile)
    shift_min = int(0 - arg_min)
    azimuthal_profile = np.roll(azimuthal_profile, shift_min)

    # Find extents with the threshold
    left_theta_i = np.searchsorted(azimuthal_profile, threshold)
    right_theta_i = len(theta) - (np.searchsorted(azimuthal_profile[::-1], threshold))

    left_theta = theta[left_theta_i]
    right_theta = theta[right_theta_i]

    extent = right_theta - left_theta
    return extent

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
