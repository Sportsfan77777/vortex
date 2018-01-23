"""
convert fits files into pickle files
"""

from astropy.io import fits
import pickle


def new_argument_parser(description = "Plot azimuthal density profiles in two by two grid."):
    parser = argparse.ArgumentParser()

    # File Selection
    parser.add_argument('fn',
                         help = 'name of imaged system')


    # Extra Properties
    parser.add_argument('-i', dest = "inclination", type = float, default = None,
                         help = 'inclination angle (default: None)')
    parser.add_argument('--pa', dest = "position_angle", type = float, default = None,
                         help = 'position angle (default: None)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()



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

	# Store in Pickle
	pickle.dump(intensity, fn1)
	pickle.dump(deprojected_intensity, fn2)
	pickle.dump(header, fn3)


### Store! ###

store()
