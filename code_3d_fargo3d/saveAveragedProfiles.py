import sys, os, subprocess
import pickle, glob
from multiprocessing import Pool
import argparse

import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot

from pylab import rcParams
from pylab import fromfile

import util
import utilVorticity
import azimuthal as az
#from readTitle import readTitle

from advanced import Parameters
#from reader_mpiio import Fields

from colormaps import cmaps
for key in cmaps:
    plot.register_cmap(name = key, cmap = cmaps[key])


def new_argument_parser(description = "Plot gas density maps."):
    parser = argparse.ArgumentParser()

    # Frame Selection
    parser.add_argument('frames', type = int, nargs = '+',
                         help = 'select single frame or range(start, end, rate). error if nargs != 1 or 3')
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')

    # Files
    parser.add_argument('--dir', dest = "save_directory", default = "averagedProfiles",
                         help = 'save directory (default: averagedProfiles)')

    # File Selection
    parser.add_argument('--gas', dest = "gas", action = 'store_true', default = False,
                         help = 'save gas-related (default: do not)')
    parser.add_argument('--dust', dest = "dust", action = 'store_true', default = False,
                         help = 'save dust-related (default: do not)')
    parser.add_argument('--vz', dest = "vz", action = 'store_true', default = False,
                         help = 'save vz-related only (default: do not)')

    # Plot Parameters (contours)

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

### Get Fargo Parameters ###
p = Parameters()

num_rad = p.ny; num_theta = p.nx; num_z = p.nz
r_min = p.ymin; r_max = p.ymax
z_min = p.zmin; z_max = p.zmax

spacing = p.spacing

surface_density_zero = p.sigma0

dt = p.ninterm * p.dt

fargo_par = util.get_pickled_parameters()
jupiter_mass = 1e-3
planet_mass = fargo_par["PlanetMass"] / jupiter_mass
accretion = fargo_par["Accretion"]
taper_time = p.masstaper

scale_height = p.aspectratio
viscosity = 1e-7 #p.nu

# Frames
frame_range = util.get_frame_range(args.frames)

# Number of Cores 
num_cores = args.num_cores

# Files
save_directory = args.save_directory
if not os.path.isdir(save_directory):
    os.mkdir(save_directory) # make save directory if it does not already exist

if spacing is "Arithmetic":
    rad = np.linspace(r_min, r_max, num_rad)
else:
    rad = np.logspace(np.log10(r_min), np.log10(r_max), num_rad)
theta = np.linspace(0, 2 * np.pi, num_theta)
z_angles = np.linspace(z_min, z_max, num_z)

###############################################################################

###### GAS ######

def save_density(frame):
    # Read data
    density = fromfile("gasdens%d.dat" % frame).reshape(num_z, num_rad, num_theta)

    # Density
    averagedDensity = np.average(density, axis = 2)
    pickle.dump(averagedDensity, open("%s/averagedDensity%04d.p" % (save_directory, frame), 'wb'))

    # Midplane density
    midplane_density = density[num_z / 2, :, :]
    pickle.dump(midplane_density, open("%s/midplaneDensity%04d.p" % (save_directory, frame), 'wb'))

    # Surface density
    dz = z_angles[1] - z_angles[0]
    surface_density = np.sum(density[:, :, :], axis = 0) * dz
    averagedSurfaceDensity = np.average(surface_density, axis = 1)

    pickle.dump(surface_density, open("%s/surfaceDensity%04d.p" % (save_directory, frame), 'wb'))
    pickle.dump(averagedSurfaceDensity, open("%s/averagedSurfaceDensity%04d.p" % (save_directory, frame), 'wb'))

def save_velocity(frame):
    # Read data
    vz = (fromfile("gasvz%d.dat" % frame).reshape(num_z, num_rad, num_theta))
    vrad = (fromfile("gasvy%d.dat" % frame).reshape(num_z, num_rad, num_theta))
    vtheta = (fromfile("gasvx%d.dat" % frame).reshape(num_z, num_rad, num_theta))

    # Average
    averaged_vz = np.average(vz, axis = 2)
    averaged_vrad = np.average(vrad, axis = 2)
    averaged_vtheta = np.average(vtheta, axis = 2)

    pickle.dump(averaged_vz, open("%s/averaged-vz%04d.p" % (save_directory, frame), 'wb'))
    pickle.dump(averaged_vrad, open("%s/averaged-vy%04d.p" % (save_directory, frame), 'wb'))
    pickle.dump(averaged_vtheta, open("%s/averaged-vx%04d.p" % (save_directory, frame), 'wb'))

    # Midplane
    midplane_vz = vz[num_z / 2, :, :]
    midplane_vrad = vrad[num_z / 2, :, :]
    midplane_vtheta = vtheta[num_z / 2, :, :]

    pickle.dump(midplane_vz, open("%s/midplane-vz%04d.p" % (save_directory, frame), 'wb'))
    pickle.dump(midplane_vrad, open("%s/midplane-vy%04d.p" % (save_directory, frame), 'wb'))
    pickle.dump(midplane_vtheta, open("%s/midplane-vx%04d.p" % (save_directory, frame), 'wb'))

def save_vz(frame):
    # Read data
    vz = (fromfile("gasvz%d.dat" % frame).reshape(num_z, num_rad, num_theta))

    # Average
    averaged_vz = np.average(vz, axis = 2)
    pickle.dump(averaged_vz, open("%s/averaged-vz%04d.p" % (save_directory, frame), 'wb'))

    # Midplane
    midplane_vz = vz[num_z / 2, :, :]
    pickle.dump(midplane_vz, open("%s/midplane-vz%04d.p" % (save_directory, frame), 'wb'))

def save_energy(frame):
    # Read data
    energy = fromfile("gasenergy%d.dat" % frame).reshape(num_z, num_rad, num_theta)

    # Density
    averagedEnergy = np.average(energy, axis = 2)
    pickle.dump(averagedEnergy, open("%s/averagedEnergy%04d.p" % (save_directory, frame), 'wb'))

    # Midplane density
    midplane_energy = energy[num_z / 2, :, :]
    pickle.dump(midplane_energy, open("%s/midplaneEnergy%04d.p" % (save_directory, frame), 'wb'))

    # Surface density
    dz = z_angles[1] - z_angles[0]
    surface_energy = np.sum(energy[:, :, :], axis = 0) * dz
    averagedSurfaceEnergy = np.average(energy, axis = 1)

    pickle.dump(surface_energy, open("%s/surfaceEnergy%04d.p" % (save_directory, frame), 'wb'))
    pickle.dump(averagedSurfaceEnergy, open("%s/averagedSurfaceEnergy%04d.p" % (save_directory, frame), 'wb'))


###############################################################################

###### DUST ######

def save_dust_density(frame):
    # Read data
    density = fromfile("dust1dens%d.dat" % frame).reshape(num_z, num_rad, num_theta)

    # Density
    averagedDensity = np.average(density, axis = 2)
    pickle.dump(averagedDensity, open("%s/averagedDust1Density%04d.p" % (save_directory, frame), 'wb'))

    # Midplane density
    midplane_density = density[num_z / 2, :, :]
    pickle.dump(midplane_density, open("%s/midplaneDust1Density%04d.p" % (save_directory, frame), 'wb'))

    # Surface density
    dz = z_angles[1] - z_angles[0]
    surface_density = np.sum(density[:, :, :], axis = 0) * dz
    averagedSurfaceDensity = np.average(surface_density, axis = 1)

    pickle.dump(surface_density, open("%s/surfaceDust1Density%04d.p" % (save_directory, frame), 'wb'))
    pickle.dump(averagedSurfaceDensity, open("%s/averagedDust1SurfaceDensity%04d.p" % (save_directory, frame), 'wb'))

###############################################################################

def save_files(frame):
    print frame

    # Gas
    if args.gas:
        save_density(frame)
        save_velocity(frame)
        save_energy(frame)

    if args.vz:
        save_vz(frame)

    # Dust
    if args.dust:
        save_dust_density(frame)

# Iterate through frames

if len(frame_range) == 1:
    save_files(frame_range[0])
else:
    if num_cores > 1:
        p = Pool(num_cores) # default number of processes is multiprocessing.cpu_count()
        p.map(save_files, frame_range)
        p.terminate()
    else:
        for frame in frame_range:
            save_files(frame)

