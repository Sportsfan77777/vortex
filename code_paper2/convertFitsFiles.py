"""
convert fits files into pickle files
"""

import argparse

from astropy.io import fits
import pickle

import numpy as np
from scipy.interpolate import RectBivariateSpline as spline


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # File Selection
    parser.add_argument('fn',
                         help = 'filename of imaged system')

    # Save Info
    parser.add_argument('--id', dest = "id_number", type = int, default = None,
                         help = 'save number (default: None)')
    parser.add_argument('--name', dest = "name", default = None,
                         help = 'name of system (default: stem of fn)')

    # Extra Properties
    parser.add_argument('-i', dest = "inclination", type = float, default = None,
                         help = 'inclination angle (default: None)')
    parser.add_argument('--pa', dest = "position_angle", type = float, default = None,
                         help = 'position angle (default: None)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

if args.name is None:
    name = args.fn[:-5] # remove .fits ending
else:
    name = args.name
basename = "fits%03d_%s.p" # savename

###############################################################################

### Helper Function ###

def deproject_image(incl, pa, image):
    npix=len(image)
    imx = 1.0 * np.arange(npix); imy = 1.0 * np.arange(npix)
    imx -= np.median(imx); imy -= np.median(imy)

    pa_rad = pa * (np.pi / 180.0)
    incl_rad = incl * (np.pi / 180.0)

    xx,yy   = np.meshgrid(imx, imy)
    xxr      = xx*np.cos(pa_rad) - yy*np.sin(pa_rad)
    yyr      = xx*np.sin(pa_rad) + yy*np.cos(pa_rad)
    spline1   = spline(imx, imy, image.T, kx=3, ky=3)
    dummy_inu = np.zeros([npix, npix], dtype=float)
    for ix in range(npix):
        for iy in range(npix):
            dummy_inu[iy,ix] = spline1(xxr[iy,ix], yyr[iy,ix])
    xxr      = (xx*np.cos(pa_rad) + yy*np.sin(pa_rad)) * np.cos(incl_rad)
    yyr      = -xx*np.sin(pa_rad) + yy*np.cos(pa_rad)
    sp       = spline(imx, imy, dummy_inu.T, kx=3, ky=3)
    inuDP    = np.zeros([npix, npix], dtype=float)
    for ix in range(npix):
        for iy in range(npix):
            inuDP[iy,ix] = sp(xxr[iy,ix], yyr[iy,ix])
    return(inuDP)

###############################################################################

def clean(header):
    """ Remove history, save all of the normal keys in a dictionary, and add fargo_par parameters """
    new_header = {}

    # Delete history and get keys
    del header['HISTORY']
    keys = header.keys()

    # Create Dictionary
    for key in keys:
        new_header[key] = header[key]

    # Add extra 'fargo_par' parameters
    px_scale = new_header['cdelt2'] * 3600
    num_x = new_header['naxis1']; num_y = new_header['naxis2']
    width = num_x * px_scale; height = num_y * px_scale

    max_r = np.sqrt(np.power(np.max(width / 2.0), 2) + np.power(np.max(height / 2.0), 2))

    new_header["rad"] = np.linspace(0, max_r, num_x)
    new_header["theta"] = np.linspace(0, 2 * np.pi, num_y)
    new_header["AspectRatio"] = 0.03
    new_header["Sigma0"] = 1

    return new_header

### Store Things ###

def store():
    # Read Fits File
    fits_file = fits.open(args.fn)[0]
    intensity = np.array(fits_file.data[0, 0, :, :]); header = clean(fits_file.header)

    # Extra Properties
    header["name"] = name

    if args.inclination is not None:
        header["inc"] = args.inclination
    if args.position_angle is not None:
        header["pa"] = args.position_angle

    # Deproject
    deprojected_intensity = deproject_image(header["inc"], header["pa"], intensity)

    # Store in Pickle
    fn1 = basename % (args.id_number, name)
    pickle.dump(intensity, open(fn1, "wb"))

    fn2 = "deprojected_%s" % fn1
    pickle.dump(deprojected_intensity, open(fn2, "wb"))

    fn3 = "params_%s" % fn1
    pickle.dump(header, open(fn3, "wb"))


### Store! ###

store()
