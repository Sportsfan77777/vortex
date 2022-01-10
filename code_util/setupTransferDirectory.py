"""
set up transfer directory to move files from tw2 to tiara
"""

import sys, os, shutil
import glob, pickle

cwd = os.getcwd()
cwd_base = cwd.split("/")[-1]

if not os.path.isdir(cwd_base):
    os.mkdir(cwd_base) # make transfer directory if it does not already exist

###############################
####### Copy Util Files #######
###############################

files_to_copy = ["twcc.sh", "jupiter.cfg", "planet0.dat", "bigplanet0.dat", "domain_x.dat", "domain_y.dat", "domain_z.dat", "tqwk0.dat", "this.out", "fargo3d"]

par_files = glob.glob("*.par")
for par_file in par_files:
    files_to_copy.append(par_file)

for fn in files_to_copy:
    shutil.copy2(fn, "%s/%s" % (cwd_base, fn))

###############################
###### Move Output Files ######
###############################

gas_files = glob.glob("gas*.dat")
#gas_density_files = glob.glob("gasdens*.dat")
#gas_vx_files = glob.glob("gasvx*.dat")
#gas_vy_files = glob.glob("gasvy*.dat")
#gas_vz_files = glob.glob("gasvz*.dat")
#gas_energy_files = glob.glob("gasenergy*.dat")

dust_files = glob.glob("dust*.dat")
#dust_density_files = glob.glob("dust1dens*.dat")
#dust_vx_files = glob.glob("dust1vx*.dat")
#dust_vy_files = glob.glob("dust1vy*.dat")
#dust_vz_files = glob.glob("dust1vz*.dat")
#dust_energy_files = glob.glob("dust1energy*.dat")

summary_files = glob.glob("summary*.dat")

for fn in gas_files:
    shutil.move(fn, "%s/%s" % (cwd_base, fn))

for fn in dust_files:
    shutil.move(fn, "%s/%s" % (cwd_base, fn))

for fn in summary_files:
    shutil.move(fn, "%s/%s" % (cwd_base, fn))
