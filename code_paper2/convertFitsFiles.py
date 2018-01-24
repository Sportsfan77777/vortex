"""
convert fits files into pickle files
"""

from astropy.io import fits
import pickle

import numpy as np
from scipy.interpolate import RectBivariateSpline as spline


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # File Selection
    parser.add_argument('fn',
                         help = 'name of imaged system')

    # Save Number
    parser.add_argument('--id', dest = "id_number", type = int, default = None,
                         help = 'save number (default: None)')

    # Extra Properties
    parser.add_argument('-i', dest = "inclination", type = float, default = None,
                         help = 'inclination angle (default: None)')
    parser.add_argument('--pa', dest = "position_angle", type = float, default = None,
                         help = 'position angle (default: None)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

name = args.fn[:-5] # remove .fits ending
basename = "fits%03d_%s.p" # savename

###############################################################################

### Helper Function ###

def deproject_image(incl, pa, image):
    npix=len(image)
    imx = 1.0 * np.arange(npix); imy = 1.0 * np.arange(npix)
    imx -= np.median(imx); imy -= np.median(imy)

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

### Store Things ###

def store():
	# Read Fits File
	fits_file = fits.open(args.fn)[0]
	intensity = fits_file.data[0, 0, :, :]; header = fits_file.header

	# Extra Properties
	if args.inclination is not None:
		header["inc"] = args.inclination
	if args.position_angle is not None:
		header["pa"] = args.position_angle

	# Deproject
	deprojected_intensity = deproject_image(header["inc"], header["pa"], intensity)

	# Store in Pickle
	fn1 = basename % (id_number, name)
	pickle.dump(intensity, fn1)

	fn2 = "deproject_%s" % fn1
	pickle.dump(deprojected_intensity, fn2)

	fn3 = "params_%s" % fn1
	pickle.dump(header, fn3)


### Store! ###

store()
