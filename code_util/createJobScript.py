"""
makes a new job script

"""

import argparse


def new_argument_parser(description = "Make a new job script."):
    parser = argparse.ArgumentParser()

    # File
    parser.add_argument("fn",
                         help = 'job file name (.sh appended to the end) that must be included, error otherwise')

    # Basic Parameters
    parser.add_argument('-c', dest = "num_cores", type = int, default = 1,
                         help = 'number of cores (default: 1)')
    parser.add_argument('-p', dest = "ptile", type = int, default = None,
                         help = 'number of cores needed on each computer (default: num_cores)')

    parser.add_argument('--err', dest = "err_name", default = "err_%I",
                         help = 'job error file name (default: err_%I)')
    parser.add_argument('--out', dest = "out_name", default = "out_%I",
                         help = 'job output file name (default: out_%I)')

    parser.add_argument('-q', dest = "queue", default = "medium",
                         help = 'queue (default: medium)')
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
if (args.ptile is None) or (args.ptile > args.num_cores):
    args.ptile = args.num_cores

###############################################################################

### Write File ###

with open(args.fn, 'w') as f:
    f.write("#!/bin/bash\n")

    ### Basic Parameters ###
    f.write("#BSUB -n %d\n" % args.num_cores)
    f.write("#BSUB -e %s\n" % args.err_name)
    f.write("#BSUB -o %s\n" % args.out_name)
    f.write('#BSUB -q "%s"\n' % args.queue)
    f.write("#BSUB -u mhammer\n")
    f.write("#BSUB -J %s\n" % args.name)
    if args.gpu:
        f.write("#BSUB -R gpu\n")
    f.write('#BSUB -R "span[ptile=%d]"\n' % args.ptile)

    # Line Break #
    f.write("\n")

    ### Modules ###
    if args.python:
        f.write("module load python/2.7.3\n")
    if args.fftw:
        f.write("module load fftw/2.1.5\n")
    if args.openmpi:
        f.write("module load openmpi\n")

    # Line Break
    f.write("\n")

    ### Job ###
    if args.mpirun:
        f.write("mpirun -np %d " % args.num_cores)
    f.write("%s " % args.job)
    f.write("> %s\n" % args.output)

    # Line Break
    f.write("\n")



