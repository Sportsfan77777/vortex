"""
makes a new job script

"""

import os
import argparse

def new_argument_parser(description = "Make a new job script."):
    parser = argparse.ArgumentParser()

    # File
    parser.add_argument("fn",
                         help = 'job file name (.sh appended to the end) that must be included, error otherwise')

    # Basic Parameters
    parser.add_argument('-n', dest = "num_nodes", type = int, default = 1,
                         help = 'number of nodes (default: 1)')
    parser.add_argument('-c', dest = "num_cpus", type = int, default = 1,
                         help = 'number of cpus per node (default: 1)')
    parser.add_argument('-t', dest = "num_hours", type = int, default = 1,
                         help = 'number of walltime hours (default: 1)')

    parser.add_argument('-m', dest = "mem", type = int, default = None,
                         help = 'memory in gb (default: 6 per cpu)')

    parser.add_argument('-q', dest = "queue", default = "standard",
                         help = 'queue (default: standard)')
    parser.add_argument('--name', dest = "name", default = None,
                         help = 'queue (default: fn)')

    parser.add_argument('--gpu', dest = "gpu", action = 'store_true', default = False,
                         help = 'request gpu resource (default: no gpus)')

    # Modules
    parser.add_argument('--python_off', dest = "python", action = 'store_false', default = True,
                         help = 'include python module (default: include)')
    parser.add_argument('--fftw_off', dest = "fftw", action = 'store_false', default = True,
                         help = 'include fftw module (default: include)')
    parser.add_argument('--openmpi_off', dest = "openmpi", action = 'store_false', default = True,
                         help = 'include openmpi module (default: include)')

    # Job
    parser.add_argument('--mpi', dest = "mpirun", action = 'store_true', default = False,
                         help = 'use mpirun (default: do not use mpirun)')
    parser.add_argument('-j', dest = "job", default = "",
                         help = 'job command (default: empty string)')
    parser.add_argument('-o', dest = "output", default = None,
                         help = 'output file (.out appended to the end) (default: name)')

    return parser

###############################################################################

### Parse Arguments ###
args = new_argument_parser().parse_args()

# Names
if args.name is None:
    args.name = args.fn
if args.output is None:
    args.output = args.name

args.fn = "%s.sh" % args.fn
args.output = "%s.out" % args.output

# Cores
args.num_cores = args.num_nodes * args.num_cpus

# Memory
if args.mem is None:
    if args.gpu:
        args.mem = 8 * args.num_cpus
    else:
        args.mem = 6 * args.num_cpus

# GPUs
args.gpu_string = ""
if args.gpu:
    args.gpu_string = ":np100s=1"

###############################################################################

### Write File ###

with open(args.fn, 'w') as f:
    f.write("#!/bin/bash\n")

    ### Basic Parameters ###
    f.write("#PBS -N %s\n" % args.name)
    f.write("#PBS -W group_list=kkratter\n")
    f.write('#PBS -q "%s"\n' % args.queue)
    f.write("#PBS -m bea\n")
    f.write("#PBS -l select=%d:ncpus=%d:mem=%dgb%s\n" % (args.num_nodes, args.num_cpus, args.mem, args.gpu_string))
    f.write("#PBS -l walltime=%d:00:00\n" % args.num_hours)
    f.write("#PBS -l cput=%d:00:00\n" % (args.num_cores * args.num_hours))

    # Line Break #
    f.write("\n")

    ### Modules ###
    if args.python:
        f.write("module load python/2.7.3\n")
    if args.fftw:
        f.write("module load fftw/2.1.5\n")
    if args.openmpi:
        f.write("module load openmpi\n")
    if args.gpu:
        f.write("module unload gcc/6.1.0\n")
        f.write("module load mpich/ge/gcc/64/3.2.1\n")
        f.write("module load cuda80/toolkit/8.0.61\n")

    # Line Break #
    f.write("\n")

    ### Set Directory
    f.write("cd %s" % os.getcwd())

    # Line Break
    f.write("\n")

    ### Job ###
    if args.mpirun:
        f.write("mpirun -np %d " % args.num_cores)
    f.write("%s " % args.job)
    f.write("> %s\n" % args.output)

    # Line Break
    f.write("\n")



